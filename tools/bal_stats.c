#include "linearizer.h"
#include "solver.h"

#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <metis.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

enum {
    POSE_DIM = 6,
    VECTOR_DIM = 3,
    CAMERA_KEYS = 2,
    RESIDUAL_DIM = 2,
    BLOCKS_PER_OBS = 6,
    MERGED_BLOCKS_PER_OBS = 3,
    CAMERA_DIM = POSE_DIM + VECTOR_DIM,
};

typedef struct {
    i32 cameras;
    i32 points;
    i32 observations;
    i32* camera_indices;
    i32* point_indices;
} bal_problem;

typedef struct {
    const char* path;
    int64_t file_bytes;
    i32 cameras;
    i32 points;
    i32 observations;
    i32 residuals;
    i32 keys;
    i32 dim;
    i32 hessian_blocks;
    i32 hessian_nnz;
    i32 factor_nnz;
    i32 fill_nnz;
    i32 seed;
    u64 factor_madds;
    u64 point_madds;
    u64 camera_madds;
    i32 schur_hessian_nnz;
    i32 schur_factor_nnz;
} bal_stats;

typedef struct {
    bal_stats baseline;
    bal_stats point_first;
} schur_stats;

typedef enum {
    REPORT_METADATA,
    REPORT_SCHUR,
} report_mode;

typedef enum {
    PARSE_OK,
    PARSE_HELP,
    PARSE_ERROR,
} parse_status;

typedef int (*stats_calc)(const char*, i32*, bal_stats*);

static void* heap_alloc(size bytes, void* context) {
    (void)context;

    void* data = malloc((usize)bytes);
    if (data == NULL) {
        fprintf(stderr, "allocation failed for %td bytes\n", bytes);
        exit(EXIT_FAILURE);
    }

    return data;
}

static void heap_free(void* data, size bytes, void* context) {
    (void)bytes;
    (void)context;

    free(data);
}

static void bal_free(bal_problem problem) {
    free(problem.camera_indices);
    free(problem.point_indices);
}

static int bal_read(const char* path, bal_problem* problem) {
    FILE* file = fopen(path, "r");
    if (file == NULL) {
        perror(path);
        return -1;
    }

    // Only observations affect the Hessian sparsity pattern.
    int read = fscanf(file, "%d %d %d", &problem->cameras, &problem->points,
                      &problem->observations);
    if (read != 3 || problem->cameras <= 0 || problem->points <= 0 ||
        problem->observations <= 0) {
        fprintf(stderr, "%s: invalid BAL header\n", path);
        fclose(file);
        return -1;
    }

    problem->camera_indices = malloc((usize)problem->observations * sizeof(i32));
    problem->point_indices = malloc((usize)problem->observations * sizeof(i32));
    if (problem->camera_indices == NULL || problem->point_indices == NULL) {
        fprintf(stderr, "%s: allocation failed\n", path);
        fclose(file);
        bal_free(*problem);
        return -1;
    }

    for (i32 i = 0; i < problem->observations; ++i) {
        double x;
        double y;
        i32 camera;
        i32 point;

        read = fscanf(file, "%d %d %lf %lf", &camera, &point, &x, &y);
        if (read != 4 || camera < 0 || camera >= problem->cameras || point < 0 ||
            point >= problem->points) {
            fprintf(stderr, "%s: invalid observation %d\n", path, i);
            fclose(file);
            bal_free(*problem);
            return -1;
        }

        problem->camera_indices[i] = camera;
        problem->point_indices[i] = point;
    }

    fclose(file);
    return 0;
}

static int64_t file_size(const char* path) {
    struct stat info;
    if (stat(path, &info) != 0) {
        return -1;
    }

    return (int64_t)info.st_size;
}

static void fill_pairs(const bal_problem* problem, i32* rows, i32* cols) {
    for (i32 i = 0; i < problem->observations; ++i) {
        i32 camera = problem->camera_indices[i];
        i32 pose = CAMERA_KEYS * camera;
        i32 intrinsics = pose + 1;
        i32 point = CAMERA_KEYS * problem->cameras + problem->point_indices[i];
        i32 start = BLOCKS_PER_OBS * i;

        rows[start + 0] = pose;
        rows[start + 1] = intrinsics;
        rows[start + 2] = point;
        rows[start + 3] = intrinsics;
        rows[start + 4] = point;
        rows[start + 5] = point;

        cols[start + 0] = pose;
        cols[start + 1] = pose;
        cols[start + 2] = pose;
        cols[start + 3] = intrinsics;
        cols[start + 4] = intrinsics;
        cols[start + 5] = point;
    }
}

static void fill_key_sizes(const bal_problem* problem, i32* sizes) {
    for (i32 i = 0; i < problem->cameras; ++i) {
        sizes[CAMERA_KEYS * i] = POSE_DIM;
        sizes[CAMERA_KEYS * i + 1] = VECTOR_DIM;
    }

    i32 point_start = CAMERA_KEYS * problem->cameras;
    for (i32 i = 0; i < problem->points; ++i) {
        sizes[point_start + i] = VECTOR_DIM;
    }
}

static u64 count_madds_range(sym_csc_mat factor, i32 begin, i32 end) {
    u64 madds = 0;
    for (i32 i = begin; i < end; ++i) {
        u64 count = (u64)(factor.col_starts[i + 1] - factor.col_starts[i]);
        madds += count * (count - 1) / 2;
    }

    return madds;
}

static u64 count_madds(sym_csc_mat factor) {
    return count_madds_range(factor, 0, factor.ncols);
}

static i32 count_factor_nnz(sym_csc_mat factor, i32 begin, i32 end) {
    i32 diagonal_nnz = end - begin;
    return diagonal_nnz + factor.col_starts[end] - factor.col_starts[begin];
}

static int calc_stats(const char* path, i32* options, bal_stats* stats) {
    bal_problem problem = {0};
    if (bal_read(path, &problem) != 0) {
        return -1;
    }

    sym_allocator allocator = {
        .malloc = heap_alloc,
        .free = heap_free,
        .ctx = NULL,
    };
    sym_allocator* alloc = &allocator;

    i32 blocks = BLOCKS_PER_OBS * problem.observations;
    i32 keys = CAMERA_KEYS * problem.cameras + problem.points;
    i32* rows = heap_alloc((size)blocks * sizeof(i32), NULL);
    i32* cols = heap_alloc((size)blocks * sizeof(i32), NULL);
    i32* block_indices = heap_alloc((size)blocks * sizeof(i32), NULL);
    fill_pairs(&problem, rows, cols);

    // Match balDemo's weighted block-level METIS ordering.
    sym_csc_mat block_hessian =
        sym_csc_from_pairs(rows, cols, blocks, keys, keys, block_indices, alloc);
    heap_free(rows, (size)blocks * sizeof(i32), NULL);
    heap_free(cols, (size)blocks * sizeof(i32), NULL);

    i32* key_sizes = heap_alloc((size)keys * sizeof(i32), NULL);
    i32* key_perm = heap_alloc((size)keys * sizeof(i32), NULL);
    fill_key_sizes(&problem, key_sizes);
    sym_get_metis_tri_perm(block_hessian, key_sizes, options, key_perm, alloc);

    // Expand blocks to the scalar Hessian, then run symbolic LDL analysis.
    sym_linearization linearization;
    sym_linearizer linearizer = sym_linearizer_new(
        block_hessian, block_indices, blocks, key_sizes, keys, key_perm, &linearization, alloc);
    i32* transpose_perm = heap_alloc((size)linearization.Hl.nnz * sizeof(i32), NULL);
    sym_csc_mat hessian_t = sym_transpose_csc(linearization.Hl, transpose_perm, alloc);
    sym_chol_factorization factorization = {0};
    sym_chol_solver solver =
        sym_new_chol_solver(hessian_t, &factorization, false, alloc);

    i32 hessian_offdiag = linearization.Hl.nnz - linearization.Hl.nrows;
    stats->path = path;
    stats->file_bytes = file_size(path);
    stats->cameras = problem.cameras;
    stats->points = problem.points;
    stats->observations = problem.observations;
    stats->residuals = RESIDUAL_DIM * problem.observations;
    stats->keys = keys;
    stats->dim = linearization.Hl.nrows;
    stats->hessian_blocks = block_hessian.nnz;
    stats->hessian_nnz = linearization.Hl.nnz;
    stats->factor_nnz = solver.L_nnz + solver.dim;
    stats->fill_nnz = solver.L_nnz - hessian_offdiag;
    stats->factor_madds = count_madds(factorization.L);

    sym_chol_solver_free(solver, alloc);
    sym_chol_factorization_free(factorization, alloc);
    sym_csc_mat_free(hessian_t, alloc);
    heap_free(transpose_perm, (size)linearization.Hl.nnz * sizeof(i32), NULL);
    sym_linearizer_free(linearizer, alloc);
    sym_linearization_free(linearization, alloc);
    sym_csc_mat_free(block_hessian, alloc);
    heap_free(key_sizes, (size)keys * sizeof(i32), NULL);
    heap_free(key_perm, (size)keys * sizeof(i32), NULL);
    heap_free(block_indices, (size)blocks * sizeof(i32), NULL);
    bal_free(problem);

    return 0;
}

static sym_csc_mat build_camera_graph(
    const bal_problem* problem, sym_allocator* alloc
) {
    i32* counts = heap_alloc((size)problem->points * sizeof(i32), NULL);
    for (i32 i = 0; i < problem->points; ++i) {
        counts[i] = 0;
    }
    for (i32 i = 0; i < problem->observations; ++i) {
        ++counts[problem->point_indices[i]];
    }

    int64_t pair_count = problem->cameras;
    for (i32 i = 0; i < problem->points; ++i) {
        pair_count += (int64_t)counts[i] * (counts[i] - 1) / 2;
    }
    if (pair_count > INT32_MAX) {
        fprintf(stderr, "camera graph exceeds 32-bit METIS limits\n");
        exit(EXIT_FAILURE);
    }

    // Group camera observations by point before forming each camera clique.
    i32* starts = heap_alloc((size)(problem->points + 1) * sizeof(i32), NULL);
    i32* cursors = heap_alloc((size)problem->points * sizeof(i32), NULL);
    i32* cameras = heap_alloc((size)problem->observations * sizeof(i32), NULL);
    starts[0] = 0;
    for (i32 i = 0; i < problem->points; ++i) {
        starts[i + 1] = starts[i] + counts[i];
        cursors[i] = starts[i];
    }
    for (i32 i = 0; i < problem->observations; ++i) {
        i32 point = problem->point_indices[i];
        cameras[cursors[point]++] = problem->camera_indices[i];
    }

    i32* rows = heap_alloc((size)pair_count * sizeof(i32), NULL);
    i32* cols = heap_alloc((size)pair_count * sizeof(i32), NULL);
    i32 pair = 0;
    for (i32 camera = 0; camera < problem->cameras; ++camera) {
        rows[pair] = camera;
        cols[pair] = camera;
        ++pair;
    }
    for (i32 point = 0; point < problem->points; ++point) {
        for (i32 i = starts[point]; i < starts[point + 1]; ++i) {
            for (i32 j = starts[point]; j < i; ++j) {
                i32 a = cameras[i];
                i32 b = cameras[j];
                rows[pair] = a > b ? a : b;
                cols[pair] = a > b ? b : a;
                ++pair;
            }
        }
    }

    sym_csc_mat graph = sym_csc_from_pairs(
        rows, cols, pair, problem->cameras, problem->cameras, NULL, alloc);

    heap_free(rows, (size)pair_count * sizeof(i32), NULL);
    heap_free(cols, (size)pair_count * sizeof(i32), NULL);
    heap_free(cameras, (size)problem->observations * sizeof(i32), NULL);
    heap_free(cursors, (size)problem->points * sizeof(i32), NULL);
    heap_free(starts, (size)(problem->points + 1) * sizeof(i32), NULL);
    heap_free(counts, (size)problem->points * sizeof(i32), NULL);

    return graph;
}

static i32 count_schur_nnz(sym_csc_mat camera_graph) {
    i32 nnz = 0;
    for (i32 col = 0; col < camera_graph.ncols; ++col) {
        for (i32 i = camera_graph.col_starts[col]; i < camera_graph.col_starts[col + 1]; ++i) {
            i32 row = camera_graph.row_indices[i];
            if (row == col) {
                nnz += CAMERA_DIM * (CAMERA_DIM + 1) / 2;
                continue;
            }

            nnz += CAMERA_DIM * CAMERA_DIM;
        }
    }

    return nnz;
}

static void fill_merged_pairs(
    const bal_problem* problem, i32* rows, i32* cols
) {
    for (i32 i = 0; i < problem->observations; ++i) {
        i32 camera = problem->camera_indices[i];
        i32 point = problem->cameras + problem->point_indices[i];
        i32 start = MERGED_BLOCKS_PER_OBS * i;

        rows[start + 0] = camera;
        rows[start + 1] = point;
        rows[start + 2] = point;

        cols[start + 0] = camera;
        cols[start + 1] = camera;
        cols[start + 2] = point;
    }
}

static int calc_point_stats(const char* path, i32* options, bal_stats* stats) {
    bal_problem problem = {0};
    if (bal_read(path, &problem) != 0) {
        return -1;
    }

    sym_allocator allocator = {
        .malloc = heap_alloc,
        .free = heap_free,
        .ctx = NULL,
    };
    sym_allocator* alloc = &allocator;

    // Point elimination creates a clique among every point's observing cameras.
    sym_csc_mat camera_graph = build_camera_graph(&problem, alloc);
    i32* camera_sizes = heap_alloc((size)problem.cameras * sizeof(i32), NULL);
    i32* camera_perm = heap_alloc((size)problem.cameras * sizeof(i32), NULL);
    for (i32 i = 0; i < problem.cameras; ++i) {
        camera_sizes[i] = CAMERA_DIM;
    }
    sym_get_metis_tri_perm(camera_graph, camera_sizes, options, camera_perm, alloc);

    i32 blocks = MERGED_BLOCKS_PER_OBS * problem.observations;
    i32 keys = problem.cameras + problem.points;
    i32* rows = heap_alloc((size)blocks * sizeof(i32), NULL);
    i32* cols = heap_alloc((size)blocks * sizeof(i32), NULL);
    i32* block_indices = heap_alloc((size)blocks * sizeof(i32), NULL);
    fill_merged_pairs(&problem, rows, cols);
    sym_csc_mat block_hessian =
        sym_csc_from_pairs(rows, cols, blocks, keys, keys, block_indices, alloc);
    heap_free(rows, (size)blocks * sizeof(i32), NULL);
    heap_free(cols, (size)blocks * sizeof(i32), NULL);

    i32* key_sizes = heap_alloc((size)keys * sizeof(i32), NULL);
    i32* key_perm = heap_alloc((size)keys * sizeof(i32), NULL);
    for (i32 i = 0; i < problem.cameras; ++i) {
        key_sizes[i] = CAMERA_DIM;
    }
    for (i32 i = 0; i < problem.points; ++i) {
        key_sizes[problem.cameras + i] = VECTOR_DIM;
        key_perm[i] = problem.cameras + i;
    }
    for (i32 i = 0; i < problem.cameras; ++i) {
        key_perm[problem.points + i] = camera_perm[i];
    }

    // Expand the forced point-first order and analyze the complete factor.
    sym_linearization linearization;
    sym_linearizer linearizer = sym_linearizer_new(
        block_hessian, block_indices, blocks, key_sizes, keys, key_perm, &linearization, alloc);
    i32* transpose_perm = heap_alloc((size)linearization.Hl.nnz * sizeof(i32), NULL);
    sym_csc_mat hessian_t = sym_transpose_csc(linearization.Hl, transpose_perm, alloc);
    sym_chol_factorization factorization = {0};
    sym_chol_solver solver = sym_new_chol_solver(hessian_t, &factorization, false, alloc);

    i32 point_dim = VECTOR_DIM * problem.points;
    i32 hessian_offdiag = linearization.Hl.nnz - linearization.Hl.nrows;
    stats->path = path;
    stats->file_bytes = file_size(path);
    stats->cameras = problem.cameras;
    stats->points = problem.points;
    stats->observations = problem.observations;
    stats->residuals = RESIDUAL_DIM * problem.observations;
    stats->keys = keys;
    stats->dim = linearization.Hl.nrows;
    stats->hessian_blocks = block_hessian.nnz;
    stats->hessian_nnz = linearization.Hl.nnz;
    stats->factor_nnz = solver.L_nnz + solver.dim;
    stats->fill_nnz = solver.L_nnz - hessian_offdiag;
    stats->point_madds = count_madds_range(factorization.L, 0, point_dim);
    stats->camera_madds = count_madds_range(factorization.L, point_dim, solver.dim);
    stats->factor_madds = stats->point_madds + stats->camera_madds;
    stats->schur_hessian_nnz = count_schur_nnz(camera_graph);
    stats->schur_factor_nnz = count_factor_nnz(factorization.L, point_dim, solver.dim);

    sym_chol_solver_free(solver, alloc);
    sym_chol_factorization_free(factorization, alloc);
    sym_csc_mat_free(hessian_t, alloc);
    heap_free(transpose_perm, (size)linearization.Hl.nnz * sizeof(i32), NULL);
    sym_linearizer_free(linearizer, alloc);
    sym_linearization_free(linearization, alloc);
    sym_csc_mat_free(block_hessian, alloc);
    sym_csc_mat_free(camera_graph, alloc);
    heap_free(key_sizes, (size)keys * sizeof(i32), NULL);
    heap_free(key_perm, (size)keys * sizeof(i32), NULL);
    heap_free(block_indices, (size)blocks * sizeof(i32), NULL);
    heap_free(camera_sizes, (size)problem.cameras * sizeof(i32), NULL);
    heap_free(camera_perm, (size)problem.cameras * sizeof(i32), NULL);
    bal_free(problem);

    return 0;
}

static const char* file_name(const char* path) {
    const char* slash = strrchr(path, '/');
    if (slash == NULL) {
        return path;
    }

    return slash + 1;
}

static void print_header(void) {
    printf("%-29s %8s %6s %7s %8s %9s %7s %7s %8s %11s %11s %9s %7s %6s %13s\n",
           "problem", "MiB", "cams", "points", "obs", "residuals", "keys", "dim",
           "H blocks", "H nnz", "L nnz", "fill", "L/H", "seed", "factor madds");
}

static void print_stats(const bal_stats* stats) {
    const double bytes_per_mib = 1024.0 * 1024.0;
    double factor_ratio = (double)stats->factor_nnz / stats->hessian_nnz;

    printf("%-29s %8.2f %6d %7d %8d %9d %7d %7d %8d %11d %11d %9d %7.3f %6d %13" PRIu64 "\n",
           file_name(stats->path), stats->file_bytes / bytes_per_mib, stats->cameras,
           stats->points, stats->observations, stats->residuals, stats->keys, stats->dim,
           stats->hessian_blocks, stats->hessian_nnz, stats->factor_nnz, stats->fill_nnz,
           factor_ratio, stats->seed, stats->factor_madds);
}

static int compare_stats(const void* left, const void* right) {
    const bal_stats* a = left;
    const bal_stats* b = right;

    if (a->cameras != b->cameras) {
        return (a->cameras > b->cameras) - (a->cameras < b->cameras);
    }

    return strcmp(a->path, b->path);
}

static int compare_schur(const void* left, const void* right) {
    const schur_stats* a = left;
    const schur_stats* b = right;

    return compare_stats(&a->baseline, &b->baseline);
}

static void print_schur_header(void) {
    printf("%-29s %7s %7s %11s %11s %11s %11s %11s %13s %13s %13s %13s %8s\n",
           "problem", "base sd", "pf sd", "base L", "pf L", "Schur H", "Schur L",
           "Schur fill", "point madds", "camera madds", "base madds", "total madds",
           "saving");
}

static void print_schur(const schur_stats* stats) {
    const bal_stats* baseline = &stats->baseline;
    const bal_stats* point_first = &stats->point_first;
    i32 schur_fill = point_first->schur_factor_nnz - point_first->schur_hessian_nnz;
    double saving = 100.0 * (1.0 - (double)point_first->factor_madds /
                                      (double)baseline->factor_madds);

    printf("%-29s %7d %7d %11d %11d %11d %11d %11d %13" PRIu64
           " %13" PRIu64 " %13" PRIu64 " %13" PRIu64 " %7.2f%%\n",
           file_name(baseline->path), baseline->seed, point_first->seed,
           baseline->factor_nnz, point_first->factor_nnz,
           point_first->schur_hessian_nnz, point_first->schur_factor_nnz, schur_fill,
           point_first->point_madds, point_first->camera_madds,
           baseline->factor_madds, point_first->factor_madds, saving);
}

static void print_usage(const char* program) {
    printf("usage: %s [OPTIONS] BAL_FILE...\n", program);
    printf("  --nseps N       separators per level\n");
    printf("  --seeds N       try seeds 0 through N-1\n");
    printf("  --pfactor N     high-degree pruning factor\n");
    printf("  --ufactor N     separator imbalance\n");
    printf("  --niter N       refinement iterations\n");
    printf("  --rtype TYPE    one or two sided refinement\n");
    printf("  --ctype TYPE    shem or rm coarsening\n");
    printf("  --no2hop        disable two-hop matching\n");
    printf("  --no-compress   disable graph compression\n");
    printf("  --ccorder       order components separately\n");
    printf("  --schur         compare forced point-first elimination\n");
}

static int numeric_option(const char* argument) {
    if (strcmp(argument, "--nseps") == 0) {
        return METIS_OPTION_NSEPS;
    }
    if (strcmp(argument, "--pfactor") == 0) {
        return METIS_OPTION_PFACTOR;
    }
    if (strcmp(argument, "--ufactor") == 0) {
        return METIS_OPTION_UFACTOR;
    }
    if (strcmp(argument, "--niter") == 0) {
        return METIS_OPTION_NITER;
    }

    return -1;
}

static int parse_i32(const char* text, i32* value) {
    errno = 0;
    char* end = NULL;
    long parsed = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || parsed < INT32_MIN ||
        parsed > INT32_MAX) {
        return -1;
    }

    *value = (i32)parsed;
    return 0;
}

static parse_status parse_args(
    int argc, char** argv, i32* options, i32* seed_count, report_mode* mode,
    char** paths, int* path_count
) {
    for (int i = 1; i < argc; ++i) {
        const char* argument = argv[i];
        if (strcmp(argument, "--help") == 0) {
            return PARSE_HELP;
        }
        if (argument[0] != '-') {
            paths[(*path_count)++] = argv[i];
            continue;
        }

        int option = numeric_option(argument);
        if (option >= 0) {
            if (++i >= argc || parse_i32(argv[i], &options[option]) != 0) {
                fprintf(stderr, "%s requires an integer\n", argument);
                return PARSE_ERROR;
            }
            continue;
        }

        if (strcmp(argument, "--seeds") == 0) {
            if (++i >= argc || parse_i32(argv[i], seed_count) != 0 || *seed_count <= 0) {
                fprintf(stderr, "--seeds requires a positive integer\n");
                return PARSE_ERROR;
            }
            continue;
        }

        if (strcmp(argument, "--schur") == 0) {
            *mode = REPORT_SCHUR;
            continue;
        }

        if (strcmp(argument, "--rtype") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "--rtype requires one or two\n");
                return PARSE_ERROR;
            }
            if (strcmp(argv[i], "one") == 0) {
                options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP1SIDED;
                continue;
            }
            if (strcmp(argv[i], "two") == 0) {
                options[METIS_OPTION_RTYPE] = METIS_RTYPE_SEP2SIDED;
                continue;
            }

            fprintf(stderr, "--rtype requires one or two\n");
            return PARSE_ERROR;
        }

        if (strcmp(argument, "--ctype") == 0) {
            if (++i >= argc) {
                fprintf(stderr, "--ctype requires shem or rm\n");
                return PARSE_ERROR;
            }
            if (strcmp(argv[i], "shem") == 0) {
                options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
                continue;
            }
            if (strcmp(argv[i], "rm") == 0) {
                options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM;
                continue;
            }

            fprintf(stderr, "--ctype requires shem or rm\n");
            return PARSE_ERROR;
        }

        if (strcmp(argument, "--no2hop") == 0) {
            options[METIS_OPTION_NO2HOP] = 1;
            continue;
        }
        if (strcmp(argument, "--no-compress") == 0) {
            options[METIS_OPTION_COMPRESS] = 0;
            continue;
        }
        if (strcmp(argument, "--ccorder") == 0) {
            options[METIS_OPTION_CCORDER] = 1;
            continue;
        }

        fprintf(stderr, "unknown option: %s\n", argument);
        return PARSE_ERROR;
    }

    return PARSE_OK;
}

static int calc_best(
    const char* path, i32* options, i32 seed_count, stats_calc calculate,
    bal_stats* best
) {
    if (seed_count == 0) {
        int result = calculate(path, options, best);
        best->seed = options[METIS_OPTION_SEED];
        return result;
    }

    for (i32 seed = 0; seed < seed_count; ++seed) {
        options[METIS_OPTION_SEED] = seed;

        bal_stats candidate;
        if (calculate(path, options, &candidate) != 0) {
            return -1;
        }
        candidate.seed = seed;

        if (seed > 0 && candidate.factor_madds > best->factor_madds) {
            continue;
        }
        if (seed > 0 && candidate.factor_madds == best->factor_madds &&
            candidate.factor_nnz >= best->factor_nnz) {
            continue;
        }

        *best = candidate;
    }

    return 0;
}

int main(int argc, char** argv) {
    i32 options[METIS_NOPTIONS];
    if (METIS_SetDefaultOptions(options) != METIS_OK) {
        fprintf(stderr, "failed to initialize METIS options\n");
        return EXIT_FAILURE;
    }

    char** paths = calloc((usize)argc, sizeof(char*));
    if (paths == NULL) {
        fprintf(stderr, "allocation failed\n");
        return EXIT_FAILURE;
    }

    int path_count = 0;
    i32 seed_count = 0;
    report_mode mode = REPORT_METADATA;
    parse_status parse =
        parse_args(argc, argv, options, &seed_count, &mode, paths, &path_count);
    if (parse == PARSE_HELP) {
        print_usage(argv[0]);
        free(paths);
        return EXIT_SUCCESS;
    }
    if (parse == PARSE_ERROR || path_count == 0) {
        print_usage(argv[0]);
        free(paths);
        return EXIT_FAILURE;
    }

    // Analyze first so the final table can be ordered by problem size.
    int result = EXIT_SUCCESS;
    int count = 0;
    if (mode == REPORT_METADATA) {
        bal_stats* stats = calloc((usize)path_count, sizeof(bal_stats));
        if (stats == NULL) {
            fprintf(stderr, "allocation failed\n");
            free(paths);
            return EXIT_FAILURE;
        }

        for (int i = 0; i < path_count; ++i) {
            if (calc_best(paths[i], options, seed_count, calc_stats, &stats[count]) != 0) {
                result = EXIT_FAILURE;
                continue;
            }

            ++count;
        }

        qsort(stats, (usize)count, sizeof(bal_stats), compare_stats);
        print_header();
        for (int i = 0; i < count; ++i) {
            print_stats(&stats[i]);
        }
        free(stats);
    } else {
        schur_stats* stats = calloc((usize)path_count, sizeof(schur_stats));
        if (stats == NULL) {
            fprintf(stderr, "allocation failed\n");
            free(paths);
            return EXIT_FAILURE;
        }

        for (int i = 0; i < path_count; ++i) {
            schur_stats* current = &stats[count];
            if (calc_best(paths[i], options, seed_count, calc_stats,
                          &current->baseline) != 0 ||
                calc_best(paths[i], options, seed_count, calc_point_stats,
                          &current->point_first) != 0) {
                result = EXIT_FAILURE;
                continue;
            }
            if (current->baseline.hessian_nnz != current->point_first.hessian_nnz) {
                fprintf(stderr, "%s: point-first Hessian pattern mismatch\n", paths[i]);
                result = EXIT_FAILURE;
                continue;
            }

            ++count;
        }

        qsort(stats, (usize)count, sizeof(schur_stats), compare_schur);
        print_schur_header();
        for (int i = 0; i < count; ++i) {
            print_schur(&stats[i]);
        }
        free(stats);
    }

    free(paths);
    return result;
}
