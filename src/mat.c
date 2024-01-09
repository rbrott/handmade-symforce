#include "mat.h"

#include <string.h>
#include <assert.h>
#include <stdio.h>

sym_vec sym_vec_new(i32 n, sym_allocator* alloc) {
    sym_vec v = {
        .n = n,
        .data = alloc->malloc(n, alloc->ctx)
    };
    return v;
}

void sym_vec_free(sym_vec v, sym_allocator* alloc) {
    alloc->free(v.data, v.n * sizeof(f64), alloc->ctx);
}

void sym_vec_zero(sym_vec v) {
    memset(v.data, 0, v.n * sizeof(f64));
}

sym_csc_mat sym_csc_mat_new(i32 nrows, i32 ncols, i32 nnz, sym_allocator* alloc) {
    sym_csc_mat mat = {
        .nrows = nrows,
        .ncols = ncols,
        .nnz = nnz,
        .col_starts = alloc->malloc((ncols + 1) * sizeof(i32), alloc->ctx),
        .row_indices = alloc->malloc(nnz * sizeof(i32), alloc->ctx),
        .data = alloc->malloc(nnz * sizeof(f64), alloc->ctx)
    };
    return mat;
}

void sym_csc_mat_free(sym_csc_mat m, sym_allocator* alloc) {
    alloc->free(m.col_starts, (m.ncols + 1) * sizeof(i32), alloc->ctx);
    alloc->free(m.row_indices, m.nnz * sizeof(i32), alloc->ctx);
    if (m.data != NULL) {
        alloc->free(m.data, m.nnz * sizeof(f64), alloc->ctx);
    }
}

void sym_csc_mat_zero(sym_csc_mat m) {
    assert(m.data != NULL);
    memset(m.data, 0, m.nnz * sizeof(f64));
}

// TODO: some callers want an inverted perm -- consider providing both or two versions?
void sym_sort_pairs(
    i32* rows, i32* cols, i32 n, 
    i32 nrows, i32 ncols, 
    i32* perm, 
    sym_allocator* alloc
) {
    i32 max = nrows > ncols ? nrows : ncols;
    i32* buf1 = alloc->malloc(max * sizeof(i32), alloc->ctx);
    i32* buf2 = alloc->malloc(n * sizeof(i32), alloc->ctx);

    for (i32 i = 0; i < nrows; ++i) {
        buf1[i] = 0;
    }
    for (i32 i = 0; i < n; ++i) {
        buf1[rows[i]] += 1;
    }
    {
        i32 psum = 0;
        for (i32 i = 0; i < nrows; ++i) {
            i32 tmp = buf1[i];
            buf1[i] = psum;
            psum += tmp;
        }
    }
    
    for (i32 i = 0; i < n; ++i) {
        buf2[buf1[rows[i]]++] = i;
    }

    for (i32 i = 0; i < ncols; ++i) {
        buf1[i] = 0;
    }
    for (i32 i = 0; i < n; ++i) {
        buf1[cols[i]] += 1;
    }
    {
        i32 psum = 0;
        for (i32 i = 0; i < ncols; ++i) {
            i32 tmp = buf1[i];
            buf1[i] = psum;
            psum += tmp;
        }
    }

    for (i32 i = 0; i < n; ++i) {
        i32 j = buf2[i];
        perm[buf1[cols[j]]++] = j;
    }

    alloc->free(buf1, max * sizeof(i32), alloc->ctx);
    alloc->free(buf2, n * sizeof(i32), alloc->ctx);
}

sym_csc_mat sym_csc_from_pairs(
    i32* rows, i32* cols, i32 n, 
    i32 nrows, i32 ncols, 
    i32* nz_indices,
    sym_allocator* alloc
) {
    i32* perm = (i32*) alloc->malloc(n * sizeof(i32), alloc->ctx);
    sym_sort_pairs(rows, cols, n, nrows, ncols, perm, alloc);

    i32 nnz = 0;
    i32* nnz_by_col = (i32*) alloc->malloc(ncols * sizeof(i32), alloc->ctx);
    for (i32 i = 0; i < ncols; ++i) {
        nnz_by_col[i] = 0;
    }
    if (nz_indices == NULL) {
        i32 last_row = -1, last_col = -1;
        for (i32 i = 0; i < n; ++i) {
            i32 j = perm[i];
            i32 row = rows[j];
            i32 col = cols[j];
            if (row != last_row || col != last_col) {
                ++nnz;
                ++nnz_by_col[col];
                last_row = row;
                last_col = col;
            }
        }
    } else {
        i32 last_row = -1, last_col = -1;
        for (i32 i = 0; i < n; ++i) {
            i32 j = perm[i];
            i32 row = rows[j];
            i32 col = cols[j];
            if (row != last_row || col != last_col) {
                ++nnz;
                ++nnz_by_col[col];
                last_row = row;
                last_col = col;
            }
            nz_indices[j] = nnz - 1;
        }
    }

    i32* col_starts = (i32*) alloc->malloc((ncols + 1) * sizeof(i32), alloc->ctx);
    col_starts[0] = 0;
    for (i32 i = 0; i < ncols; ++i) {
        col_starts[i + 1] = col_starts[i] + nnz_by_col[i];
    }

    alloc->free(nnz_by_col, ncols * sizeof(i32), alloc->ctx);

    i32* row_indices = (i32*) alloc->malloc(nnz * sizeof(i32), alloc->ctx);
    {
        i32 nz_index = 0;
        i32 last_row = -1, last_col = -1;
        for (i32 i = 0; i < n; ++i) {
            i32 j = perm[i];
            i32 row = rows[j];
            i32 col = cols[j];
            if (row != last_row || col != last_col) {
                row_indices[nz_index++] = row;
                last_row = row;
                last_col = col;
            }
        }
    }

    alloc->free(perm, n * sizeof(i32), alloc->ctx);

    sym_csc_mat m = {
        .nrows = nrows,
        .ncols = ncols,
        .nnz = nnz,
        .col_starts = col_starts,
        .row_indices = row_indices,
        .data = NULL,
    };
    return m;
}

sym_csc_mat sym_csc_from_deduped_pairs(
    i32* rows, i32* cols, i32 n, 
    i32 nrows, i32 ncols, 
    i32* perm,
    sym_allocator* alloc
) {
    sym_sort_pairs(rows, cols, n, nrows, ncols, perm, alloc);

    i32 nnz = 0;
    i32* nnz_by_col = (i32*) alloc->malloc(ncols * sizeof(i32), alloc->ctx);
    for (i32 i = 0; i < ncols; ++i) {
        nnz_by_col[i] = 0;
    }
    for (i32 i = 0; i < n; ++i) {
        i32 j = perm[i];
        i32 row = rows[j];
        i32 col = cols[j];
        ++nnz;
        ++nnz_by_col[col];
    }

    i32* col_starts = (i32*) alloc->malloc((ncols + 1) * sizeof(i32), alloc->ctx);
    col_starts[0] = 0;
    for (i32 i = 0; i < ncols; ++i) {
        col_starts[i + 1] = col_starts[i] + nnz_by_col[i];
    }

    alloc->free(nnz_by_col, ncols * sizeof(i32), alloc->ctx);

    i32* row_indices = (i32*) alloc->malloc(nnz * sizeof(i32), alloc->ctx);
    {
        i32 nz_index = 0;
        for (i32 i = 0; i < n; ++i) {
            i32 j = perm[i];
            i32 row = rows[j];
            i32 col = cols[j];
            row_indices[nz_index++] = row;
        }
    }

    sym_csc_mat m = {
        .nrows = nrows,
        .ncols = ncols,
        .nnz = nnz,
        .col_starts = col_starts,
        .row_indices = row_indices,
        .data = NULL,
    };
    return m;
}

void sym_print_csc_mat(sym_csc_mat m, sym_allocator* alloc) {
    i32* rows = (i32*) alloc->malloc(m.nnz * sizeof(i32), alloc->ctx);
    i32* cols = (i32*) alloc->malloc(m.nnz * sizeof(i32), alloc->ctx); 
    for (i32 i = 0; i < m.ncols; ++i) {
        for (i32 j = m.col_starts[i]; j < m.col_starts[i + 1]; ++j) {
            rows[j] = i;
            cols[j] = m.row_indices[j];
        }
    }

    i32* perm = (i32*) alloc->malloc(m.nnz * sizeof(i32), alloc->ctx);
    sym_csc_mat mt = sym_csc_from_deduped_pairs(rows, cols, m.nnz, m.nrows, m.ncols, perm, alloc);
    alloc->free(rows, m.nnz * sizeof(i32), alloc->ctx);
    alloc->free(cols, m.nnz * sizeof(i32), alloc->ctx);

    f64* data = m.data == NULL ? NULL : (f64*) alloc->malloc(m.nnz * sizeof(f64), alloc->ctx);
    if (data != NULL) {
        for (i32 i = 0; i < m.nnz; ++i) {
            data[i] = m.data[perm[i]];
        }
    }
    mt.data = data;

    i32 nz_index = 0;
    for (i32 i = 0; i < m.nrows; ++i) {
        for (i32 j = 0; j < m.ncols; ++j) {
            if (nz_index < mt.col_starts[i + 1] && mt.row_indices[nz_index] == j) {
                if (mt.data == NULL) {
                    printf("%d ", perm[nz_index]);
                } else {
                    printf("%f ", mt.data[nz_index]);
                }
                ++nz_index;
            } else {
                printf("_ ");
            }
        }
        printf("\n");
        assert(nz_index == mt.col_starts[i + 1]);
    }

    alloc->free(perm, m.nnz * sizeof(i32), alloc->ctx);

    sym_csc_mat_free(mt, alloc);
}
