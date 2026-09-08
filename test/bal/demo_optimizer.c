/* ----------------------------------------------------------------------------
 * SymForce - Copyright 2022, Skydio, Inc.
 * This source code is under the Apache 2.0 license found in the LICENSE file.
 * ---------------------------------------------------------------------------- */

#include <float.h>
#include <stdbool.h>
#include <string.h>

#include "gen/snavely_reprojection_factor.h"
#include "gen/pose3_retract.h"
#include "gen/rot3_tangent.h"

#include "optimizer.h"
#include "sym_assert.h"
#include "alloc.h"
#include "arena.h"
#include "linearizer.h"
#include "solver.h"
#include "types.h"

typedef struct {
    i32* camera_indices;
    i32* point_indices;
    f64* pixels;
    f64* values;
    i32 num_cameras;
    i32 num_points;
    i32 num_observations;
} bal_problem;

typedef struct {
    bal_problem problem;

    sym_linearizer lzr;

    i32* Hlt_perm;
    sym_csc_mat Hlt;

    sym_chol_factorization fac;
    sym_chol_solver solver;
} bal_optimizer_state;

bal_problem bal_read_new(char* filepath, sym_allocator* alloc) {
    FILE* file = fopen(filepath, "r");

    bal_problem p = {};
    fscanf(file, "%d", &p.num_cameras);
    fscanf(file, "%d", &p.num_points);
    fscanf(file, "%d", &p.num_observations);

    p.camera_indices = (i32*) alloc->malloc(p.num_observations * sizeof(i32), alloc->ctx);
    p.point_indices = (i32*) alloc->malloc(p.num_observations * sizeof(i32), alloc->ctx);
    // x, y interleaved
    p.pixels = (f64*) alloc->malloc(2 * p.num_observations * sizeof(f64), alloc->ctx);

    for (i32 i = 0; i < p.num_observations; i++) {
        i32 camera, point;
        fscanf(file, "%d", &camera);
        fscanf(file, "%d", &point);

        f64 px, py;
        fscanf(file, "%lf", &px);
        fscanf(file, "%lf", &py);

        p.camera_indices[i] = camera;
        p.point_indices[i] = point;
        p.pixels[2 * i + 0] = px;
        p.pixels[2 * i + 1] = py;
    }

    f64 epsilon = 10.0 * DBL_EPSILON;

    // (cam_T_world: 7, intrinsics: 3, point: 3)
    i32 values_dim = (7 + 3) * p.num_cameras + 3 * p.num_points;
    p.values = (f64*) alloc->malloc(values_dim * sizeof(f64), alloc->ctx);
    for (i32 i = 0; i < p.num_cameras; i++) {
        f64 rx, ry, rz, tx, ty, tz, f, k1, k2;
        fscanf(file, "%lf", &rx);
        fscanf(file, "%lf", &ry);
        fscanf(file, "%lf", &rz);
        fscanf(file, "%lf", &tx);
        fscanf(file, "%lf", &ty);
        fscanf(file, "%lf", &tz);
        fscanf(file, "%lf", &f);
        fscanf(file, "%lf", &k1);
        fscanf(file, "%lf", &k2);

        f64 rot_tangent[3] = {rx, ry, rz};
        rot3_from_tangent(rot_tangent, p.values + 10 * i, epsilon);

        p.values[10 * i + 4] = tx;
        p.values[10 * i + 5] = ty;
        p.values[10 * i + 6] = tz;

        p.values[10 * i + 7] = f;
        p.values[10 * i + 8] = k1;
        p.values[10 * i + 9] = k2;
    }

    for (i32 i = 0; i < p.num_points; i++) {
        f64 x, y, z;
        fscanf(file, "%lf", &x);
        fscanf(file, "%lf", &y);
        fscanf(file, "%lf", &z);

        p.values[10 * p.num_cameras + 3 * i + 0] = x;
        p.values[10 * p.num_cameras + 3 * i + 1] = y;
        p.values[10 * p.num_cameras + 3 * i + 2] = z;
    }

    fclose(file);

    return p;
}

void bal_free(bal_problem p, sym_allocator* alloc) {
    alloc->free(p.camera_indices, p.num_observations * sizeof(i32), alloc->ctx);
    alloc->free(p.point_indices, p.num_observations * sizeof(i32), alloc->ctx);
    alloc->free(p.point_indices, 2 * p.num_observations * sizeof(f64), alloc->ctx);

    i32 values_dim = (7 + 3) * p.num_cameras + 3 * p.num_points;
    alloc->free(p.values, values_dim * sizeof(f64), alloc->ctx);
}


f64 bal_linearize(
    sym_vec state,
    sym_linearization lin,
    void* userdata
) {
    bal_optimizer_state* opt_state = (bal_optimizer_state*) userdata;

    // Run a single linearization round.
    sym_linearization_clear(lin);

    // (cam_T_world: 6, intrinsics: 3, point: 3)
    f64 fac_res[2];
    f64 fac_hessian_dense[12 * 12];
    f64 fac_rhs[12];

    f64 error = 0.0;
    for (i32 obs_index = 0; obs_index < opt_state->problem.num_observations; obs_index++) {
        i32 camera_index = opt_state->problem.camera_indices[obs_index];
        i32 pose_key = 2 * camera_index + 0;
        i32 intrinsics_key = 2 * camera_index + 1;
        i32 point_index = opt_state->problem.point_indices[obs_index];
        i32 point_key = 2 * opt_state->problem.num_cameras + point_index;

        // 100+ have problems with camera at index 72
        // 257 also has a bad camera at index 238
        // TODO: focal length seems to be a problem
        //   this could probably be fixed with a prior or different focal length parameterization
        // if (camera_index == 72 || camera_index == 238) {
        //   continue;
        // }

        f64 epsilon = 10.0 * DBL_EPSILON;
        snavely_reprojection_factor(state.data + 10 * camera_index,
                                    state.data + 10 * camera_index + 7,
                                    state.data + 10 * opt_state->problem.num_cameras + 3 * point_index,
                                    opt_state->problem.pixels + 2 * obs_index,
                                    epsilon,
                                    fac_res, NULL, fac_hessian_dense, fac_rhs);

        error += fac_res[0] * fac_res[0];
        error += fac_res[1] * fac_res[1];

        // NOTE: The block index order must match the order in the block triplets above.
        sym_linearizer_add_hessian_tri_block(
            opt_state->lzr, lin,
            6 * obs_index + 0, pose_key,
            fac_hessian_dense, 12, 0
        );
        sym_linearizer_add_hessian_rect_block(
            opt_state->lzr, lin,
            6 * obs_index + 1, intrinsics_key, pose_key,
            fac_hessian_dense, 12, 6, 0
        );
        sym_linearizer_add_hessian_rect_block(
            opt_state->lzr, lin,
            6 * obs_index + 2, point_key, pose_key,
            fac_hessian_dense, 12, 9, 0
        );

        sym_linearizer_add_hessian_tri_block(
            opt_state->lzr, lin,
            6 * obs_index + 3, intrinsics_key,
            fac_hessian_dense, 12, 6
        );
        sym_linearizer_add_hessian_rect_block(
            opt_state->lzr, lin,
            6 * obs_index + 4, point_key, intrinsics_key,
            fac_hessian_dense, 12, 9, 6
        );

        sym_linearizer_add_hessian_tri_block(
            opt_state->lzr, lin,
            6 * obs_index + 5, point_key,
            fac_hessian_dense, 12, 9
        );

        sym_linearizer_add_rhs_block(
            opt_state->lzr, lin, pose_key,
            fac_rhs, 0
        );
        sym_linearizer_add_rhs_block(
            opt_state->lzr, lin, intrinsics_key,
            fac_rhs, 6
        );
        sym_linearizer_add_rhs_block(
            opt_state->lzr, lin, point_key,
            fac_rhs, 9
        );
    }
    return 0.5 * error;
}

void bal_solve(
    sym_linearization lin,
    sym_vec delta,
    void* userdata
) {
    bal_optimizer_state* opt_state = (bal_optimizer_state*) userdata;

    for (i32 i = 0; i < lin.Hl.nnz; ++i) {
        opt_state->Hlt.data[i] = lin.Hl.data[opt_state->Hlt_perm[i]];
    }

    sym_chol_solver_factor(opt_state->solver, opt_state->Hlt, opt_state->fac);

    for (i32 i = 0; i < lin.Hl.nrows; ++i) {
        delta.data[i] = -lin.rhs.data[i];
    }

    sym_chol_solver_solve_in_place(opt_state->fac, delta);
}

void bal_retract(
    sym_vec state, sym_vec delta, sym_vec new_state, void* userdata
) {
    bal_optimizer_state* opt_state = (bal_optimizer_state*) userdata;

    memcpy(new_state.data, state.data, state.n * sizeof(f64));

    f64 epsilon = 10.0 * DBL_EPSILON;

    for (i32 i = 0; i < opt_state->problem.num_cameras; ++i) {
        i32 pose_key = 2 * i + 0;
        i32 pose_values_offset = 10 * i;
        i32 pose_rhs_offset = opt_state->lzr.key_size_scan[opt_state->lzr.key_iperm[pose_key]];
        sym_pose3_retract_in_place(new_state.data + pose_values_offset, delta.data + pose_rhs_offset, epsilon);

        i32 intrinsics_key = 2 * i + 1;
        i32 intrinsics_values_offset = 10 * i + 7;
        i32 intrinsics_rhs_offset = opt_state->lzr.key_size_scan[opt_state->lzr.key_iperm[intrinsics_key]];
        for (i32 j = 0; j < 3; ++j) {
            new_state.data[intrinsics_values_offset + j] += delta.data[intrinsics_rhs_offset + j];
        }
    }
    for (i32 i = 0; i < opt_state->problem.num_points; ++i) {
        i32 point_key = 2 * opt_state->problem.num_cameras + i;
        i32 point_values_offset = 10 * opt_state->problem.num_cameras + 3 * i;
        i32 point_rhs_offset = opt_state->lzr.key_size_scan[opt_state->lzr.key_iperm[point_key]];
        for (i32 j = 0; j < 3; ++j) {
            new_state.data[point_values_offset + j] += delta.data[point_rhs_offset + j];
        }
    }
}

int main(int argc, char** argv) {
    SYM_ASSERT(argc == 2 || argc == 3);

    bool populate_lt = argc == 3;
    if (populate_lt) {
        SYM_ASSERT(strcmp(argv[2], "--populate-lt") == 0);
    }

    f64 epsilon = 10.0 * DBL_EPSILON;

    size n = 1L << 30; // 1 GiB seems fine
    printf("Allocating %ld bytes\n", n);
    u8* buf = (u8*) malloc(n);

    sym_arena arena = {
        .beg = buf,
        .end = buf + n,
    };

    sym_allocator alloc_struct = {
        .malloc = sym_arena_malloc,
        .free = sym_arena_free,
        .ctx = &arena
    };
    sym_allocator* alloc = &alloc_struct;

    bal_problem p = bal_read_new(argv[1], alloc);

    i32 nblocks = p.num_observations * 6;
    i32 nkeys = 2 * p.num_cameras + p.num_points;

    i32* Hl_block_nz_indices;
    sym_linearization lin;
    sym_linearizer lzr;
    {
        // Compute Hessian_lower block triplets.
        i32* Hl_block_rows = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
        i32* Hl_block_cols = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);

        // Linearizer callers are responsible for choosing the order of the keys.
        // Here we have all the cameras in order (pose then intrinsics) followed by all the points.
        for (i32 obs_index = 0; obs_index < p.num_observations; ++obs_index) {
            i32 camera_index = p.camera_indices[obs_index];
            i32 pose_key = 2 * camera_index + 0;
            i32 intrinsics_key = 2 * camera_index + 1;
            i32 point_index = p.point_indices[obs_index];
            i32 point_key = 2 * p.num_cameras + point_index;

            // NOTE: The order here must be consistent with the order later on in the update calls.
            Hl_block_rows[6 * obs_index + 0] = pose_key;
            Hl_block_rows[6 * obs_index + 1] = intrinsics_key;
            Hl_block_rows[6 * obs_index + 2] = point_key;
            Hl_block_rows[6 * obs_index + 3] = intrinsics_key;
            Hl_block_rows[6 * obs_index + 4] = point_key;
            Hl_block_rows[6 * obs_index + 5] = point_key;

            Hl_block_cols[6 * obs_index + 0] = pose_key;
            Hl_block_cols[6 * obs_index + 1] = pose_key;
            Hl_block_cols[6 * obs_index + 2] = pose_key;
            Hl_block_cols[6 * obs_index + 3] = intrinsics_key;
            Hl_block_cols[6 * obs_index + 4] = intrinsics_key;
            Hl_block_cols[6 * obs_index + 5] = point_key;
        }

        // Create the linearizer, linearization.
        Hl_block_nz_indices = (i32*) alloc->malloc(nblocks * sizeof(i32), alloc->ctx);
        sym_csc_mat Hl_block = sym_csc_from_pairs(Hl_block_rows, Hl_block_cols, nblocks, nkeys, nkeys, Hl_block_nz_indices, alloc);
        alloc->free(Hl_block_rows, nblocks * sizeof(i32), alloc->ctx);
        alloc->free(Hl_block_cols, nblocks * sizeof(i32), alloc->ctx);

        // Compute key sizes.
        i32* key_sizes = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
        for (i32 i = 0; i < p.num_cameras; ++i) {
            key_sizes[2 * i + 0] = 6;
            key_sizes[2 * i + 1] = 3;
        }
        for (i32 i = 0; i < p.num_points; ++i) {
            key_sizes[2 * p.num_cameras + i] = 3;
        }

        i32* key_perm = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
        sym_get_metis_tri_perm(Hl_block, key_sizes, NULL, key_perm, alloc);

        lzr = sym_linearizer_new(
            Hl_block, Hl_block_nz_indices, nblocks,
            key_sizes, nkeys,
            key_perm,
            &lin,
            alloc
        );

        sym_csc_mat_free(Hl_block, alloc);

        alloc->free(key_perm, nkeys * sizeof(i32), alloc->ctx);
        alloc->free(key_sizes, nkeys * sizeof(i32), alloc->ctx);
    }

    i32* Hlt_perm = (i32*) alloc->malloc(lin.Hl.nnz * sizeof(i32), alloc->ctx);
    sym_csc_mat Hlt = sym_transpose_csc(lin.Hl, Hlt_perm, alloc);
    Hlt.data = (f64*) alloc->malloc(lin.Hl.nnz * sizeof(f64), alloc->ctx);

    sym_chol_factorization fac = {};
    sym_chol_solver solver = sym_new_chol_solver(Hlt, &fac, populate_lt, alloc);

    bal_optimizer_state bal_opt_state = {
        .problem = p,
        .lzr = lzr,
        .Hlt_perm = Hlt_perm,
        .Hlt = Hlt,
        .fac = fac,
        .solver = solver
    };

    i32 values_dim = (7 + 3) * p.num_cameras + 3 * p.num_points;
    sym_vec values_vec = {
        .n = values_dim,
        .data = p.values
    };
    sym_optimizer* opt = sym_optimizer_new(
        bal_linearize, bal_solve, bal_retract,
        values_vec, lin, alloc, &bal_opt_state
    );

    // this is where the optimizer begins
    sym_optimizer_status opt_status;
    while ((opt_status = sym_optimizer_step(opt)) == IN_PROGRESS);

    sym_optimizer_free(opt, alloc);

    sym_chol_solver_free(solver, alloc);
    sym_chol_factorization_free(fac, alloc);

    sym_csc_mat_free(Hlt, alloc);
    alloc->free(Hlt_perm, lin.Hl.nnz * sizeof(i32), alloc->ctx);

    sym_linearizer_free(lzr, alloc);
    sym_linearization_free(lin, alloc);

    alloc->free(Hl_block_nz_indices, nblocks * sizeof(i32), alloc->ctx);

    bal_free(p, alloc);

    printf("nalloc = %ld\n", arena.nalloc);
    printf("max_nalloc = %ld\n", arena.max_nalloc);
    SYM_ASSERT(arena.nalloc == 0);

    free(buf);
}
