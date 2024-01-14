#include "linearizer.h"
#include "alloc.h"
#include "arena.h"
#include "solver.h"
#include "sym_assert.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>


int main() {
    // use one big arena for tests
    size n = 1 << 30; // 1 GiB seems fine
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

    // test sym_sort_pairs
    printf("=== Testing sym_sort_pairs ===\n");
    {
        i32 rows[3] = {2, 1, 0};
        i32 cols[3] = {0, 1, 0};
        i32 perm[3];
        sym_sort_pairs(rows, cols, 3, 3, 3, perm, alloc);
        assert(perm[0] == 2);
        assert(perm[1] == 0);
        assert(perm[2] == 1);
    }
    {
        i32 rows[5] = {3, 2, 1, 0, 0};
        i32 cols[5] = {0, 1, 0, 1, 0};
        i32 perm[5];
        sym_sort_pairs(rows, cols, 5, 4, 4, perm, alloc);
        assert(perm[0] == 4);
        assert(perm[1] == 2);
        assert(perm[2] == 0);
        assert(perm[3] == 3);
        assert(perm[4] == 1);
    }
    assert(arena.nalloc == 0);

    // test sym_linearizer
    printf("=== Testing sym_linearizer ===\n");
    {
        i32 nkeys = 3;
        i32 key_sizes[3] = {1, 2, 3};
        i32 nnz = 6;
        i32 nblocks = 2 * nnz;
        i32 Hl_block_rows[6 * 2] = {2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 0, 0};
        i32 Hl_block_cols[6 * 2] = {2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};
        i32 Hl_block_nz_indices[6 * 2];
        sym_csc_mat Hl_block = sym_csc_from_pairs(Hl_block_rows, Hl_block_cols, nblocks, nkeys, nkeys, Hl_block_nz_indices, alloc);

        // check Hl_block_nz_indices
        for (i32 i = 0; i < nblocks; i++) {
            assert(Hl_block_nz_indices[i] == (nblocks - 1 - i) / 2);
        }

        i32 key_perm[3] = {2, 0, 1};

        sym_linearization lin;
        sym_linearizer lzr = sym_linearizer_new(
            Hl_block, Hl_block_nz_indices, nblocks,
            key_sizes, nkeys, 
            key_perm,
            &lin,
            alloc
        );

        // check Hl_block_nz_iperm
        {
            i32 exp_Hl_block_nz_iperm[6] = {3, 4, 1, 5, 2, 0};
            for (i32 i = 0; i < 6; i++) {
                assert(lzr.Hl_block_nz_iperm[i] == exp_Hl_block_nz_iperm[i]);
            }
        }

        // check Hl_key_starts_by_block_nz
        { 
            i32 exp_Hl_key_starts_by_block_nz[6] = {0, 3, 4, 0, 1, 0};
            for (i32 i = 0; i < 6; i++) {
                assert(lzr.Hl_key_starts_by_block_nz[i] == exp_Hl_key_starts_by_block_nz[i]);
            }
        }

        sym_linearization_clear(lin);
        f64 ones[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        for (i32 i = 0; i < nblocks; i++) {
            i32 row_key = Hl_block_rows[i];
            i32 col_key = Hl_block_cols[i];

            if (row_key == col_key) {
                i32 size = key_sizes[row_key];
                sym_linearizer_add_hessian_tri_block(lzr, lin, i, row_key, ones, size, 0);
            } else {
                i32 row_size = key_sizes[row_key];
                i32 col_size = key_sizes[col_key];
                sym_linearizer_add_hessian_rect_block(lzr, lin, i, row_key, col_key, ones, col_size, 1, 0);
            }
        }

        // check lin.Hl.data
        for (i32 i = 0; i < lin.Hl.nnz; i++) {
            assert(lin.Hl.data[i] == 2.0);
        }

        sym_csc_mat_free(Hl_block, alloc);

        sym_linearizer_free(lzr, alloc);
        sym_linearization_free(lin, alloc);

        assert(arena.nalloc == 0);
    }

    {
        i32 nkeys = 3;
        i32 key_sizes[3] = {1, 2, 3};
        i32 nnz = 5;
        i32 nblocks = 2 * nnz;
        // missing (0, 2)
        i32 Hl_block_rows[5 * 2] = {2, 2, 2, 2, 1, 1, 1, 1, 0, 0};
        i32 Hl_block_cols[5 * 2] = {2, 2, 1, 1, 1, 1, 0, 0, 0, 0};
        i32 Hl_block_nz_indices[5 * 2];
        sym_csc_mat Hl_block = sym_csc_from_pairs(Hl_block_rows, Hl_block_cols, nblocks, nkeys, nkeys, Hl_block_nz_indices, alloc);

        i32 key_perm[3] = {2, 0, 1};

        sym_linearization lin;
        sym_linearizer lzr = sym_linearizer_new(
            Hl_block, Hl_block_nz_indices, nblocks,
            key_sizes, nkeys, 
            key_perm,
            &lin,
            alloc
        );

        // check Hl_block_nz_iperm
        {
            i32 exp_Hl_block_nz_iperm[5] = { 2, 3, 4, 1, 0 };
            for (i32 i = 0; i < 5; i++) {
                assert(lzr.Hl_block_nz_iperm[i] == exp_Hl_block_nz_iperm[i]);
            }
        }

        sym_linearization_clear(lin);
        f64 ones[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        for (i32 i = 0; i < nblocks; i++) {
            i32 row_key = Hl_block_rows[i];
            i32 col_key = Hl_block_cols[i];

            if (row_key == col_key) {
                i32 size = key_sizes[row_key];
                sym_linearizer_add_hessian_tri_block(lzr, lin, i, row_key, ones, size, 0);
            } else {
                i32 row_size = key_sizes[row_key];
                i32 col_size = key_sizes[col_key];
                sym_linearizer_add_hessian_rect_block(lzr, lin, i, row_key, col_key, ones, col_size, 1, 0);
            }
        }

        // check lin.Hl.data
        for (i32 i = 0; i < lin.Hl.nnz; i++) {
            assert(lin.Hl.data[i] == 2.0);
        }

        sym_csc_mat_free(Hl_block, alloc);

        sym_linearizer_free(lzr, alloc);
        sym_linearization_free(lin, alloc);

        assert(arena.nalloc == 0);
    }

    // test sym_solver
    printf("=== Testing sym_solver ===\n");
    {
        // TODO: this isn't quite a proper factorization, but we'll go with it anyway
        // L = {{0, 0}, {0, 0.5}}
        // D = diag({4, 1})
        // L D L^t = {{4, 2}, {2, 2}}
        i32 L_col_starts[2] = {0, 1};
        i32 L_row_indices[1] = {1};
        f64 L_data[1] = {0.5};
        f64 D_data[2] = {4.0, 1.0};
        sym_chol_factorization fac = {
            .L = {
                .nrows = 2,
                .ncols = 2,
                .nnz = 1,
                .col_starts = L_col_starts,
                .row_indices = L_row_indices,
                .data = L_data,
            },
            .D = {
                .n = 2,
                .data = D_data,
            }
        };

        f64 x_data[2] = {4.0, 2.0};
        sym_vec x = {
            .n = 2,
            .data = x_data,
        };
        sym_chol_solver_solve_in_place(fac, x, alloc);
        assert(x.data[0] == 1.0);
        assert(x.data[1] == 0.0);

        x.data[0] = 2.0;
        x.data[1] = 0.0;
        sym_chol_solver_solve_in_place(fac, x, alloc);
        assert(x.data[0] == 1.0);
        assert(x.data[1] == -1.0);

        i32 A_col_starts[3] = {0, 1, 3};
        i32 A_row_indices[3] = {0, 0, 1};
        f64 A_data[3] = {4.0, 2.0, 2.0};
        sym_csc_mat A = {
            .nrows = 2,
            .ncols = 2,
            .nnz = 3,
            .col_starts = A_col_starts,
            .row_indices = A_row_indices,
            .data = A_data,
        };
        sym_chol_solver solver = sym_new_chol_solver(A, &fac, alloc);
        sym_chol_solver_factor(solver, A, fac);
        assert(fac.L.nnz == 1);
        assert(fac.L.col_starts[0] == 0);
        assert(fac.L.col_starts[1] == 1);
        assert(fac.L.row_indices[0] == 1);
        assert(fac.L.data[0] == 0.5);

        sym_chol_solver_free(solver, alloc);
        sym_chol_factorization_free(fac, alloc);

        assert(arena.nalloc == 0);
    }

    // test sym_transpose_csc
    printf("=== Testing sym_transpose_csc ===\n");
    {
        i32 A_col_starts[4] = {0, 1, 2, 4};
        i32 A_row_indices[4] = {1, 1, 0, 1};
        sym_csc_mat A = {
            .nrows = 2,
            .ncols = 3,
            .nnz = 4,
            .col_starts = A_col_starts,
            .row_indices = A_row_indices,
            .data = NULL,
        };
        i32 perm[4];
        sym_csc_mat At = sym_transpose_csc(A, perm, alloc);

        assert(A.nrows == At.ncols);
        assert(A.ncols == At.nrows);
        assert(A.nnz == At.nnz);

        assert(At.col_starts[0] == 0);
        assert(At.col_starts[1] == 1);
        assert(At.col_starts[2] == 4);

        assert(At.row_indices[0] == 2);
        assert(At.row_indices[1] == 0);
        assert(At.row_indices[2] == 1);
        assert(At.row_indices[3] == 2);

        f64 A_data[4] = {1.0, 2.0, 3.0, 4.0};
        f64 At_data[4];
        for (i32 i = 0; i < 4; ++i) {
            At_data[i] = A_data[perm[i]];
        }

        assert(At_data[0] == 3.0);
        assert(At_data[1] == 1.0);
        assert(At_data[2] == 2.0);
        assert(At_data[3] == 4.0);

        sym_csc_mat_free(At, alloc);

        assert(arena.nalloc == 0);
    }
}