#include "linearizer.h"
#include "alloc.h"
#include "arena.h"

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
}