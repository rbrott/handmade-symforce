#include "linearizer.h"

#include <stdio.h>
#include <assert.h>

#include <metis.h>

void sym_get_metis_tri_perm(sym_csc_mat m, i32* weights, i32* perm, sym_allocator* alloc) {
    assert(m.nrows == m.ncols);

    // METIS wants a symmetric matrix with no entries on the diagonal 
    i32 ndiag = 0;
    for (i32 i = 0; i < m.ncols; ++i) {
        if (m.col_starts[i] == m.col_starts[i + 1]) {
            continue;
        }
        ndiag += (m.row_indices[m.col_starts[i]] == i);
    }

    i32 m2_nnz = 2 * (m.nnz - ndiag);
    i32* m2_rows = (i32*) alloc->malloc(m2_nnz * sizeof(i32), alloc->ctx);
    i32* m2_cols = (i32*) alloc->malloc(m2_nnz * sizeof(i32), alloc->ctx);
    {
      i32 col = 0;
      i32 j = 0;
      for (i32 i = 0; i < m.nnz; ++i) {
        while (m.col_starts[col + 1] <= i) {
          ++col;
        }
        i32 row = m.row_indices[i];

        if (row == col) {
          continue;
        }

        m2_rows[j] = row;
        m2_cols[j] = col;

        m2_rows[j+1] = col;
        m2_cols[j+1] = row;

        j += 2;
      }
      assert(j == m2_nnz);
    }

    i32* temp_perm = (i32*) alloc->malloc(m2_nnz * sizeof(i32), alloc->ctx);
    sym_csc_mat m2 = sym_csc_from_deduped_pairs(m2_rows, m2_cols, m2_nnz, m.nrows, m.ncols, temp_perm, alloc);
    assert(m2.nnz == m2_nnz);

    alloc->free(temp_perm, m2_nnz * sizeof(i32), alloc->ctx);

    alloc->free(m2_rows, m2_nnz * sizeof(i32), alloc->ctx);
    alloc->free(m2_cols, m2_nnz * sizeof(i32), alloc->ctx);

    i32* iperm = (i32*) alloc->malloc(m.ncols * sizeof(i32), alloc->ctx); // yep, this is required
    int result = METIS_NodeND(&m.ncols, m2.col_starts, m2.row_indices, weights, NULL, perm, iperm);
    assert(result == METIS_OK);
    alloc->free(iperm, m.ncols * sizeof(i32), alloc->ctx);
    sym_csc_mat_free(m2, alloc);
}

sym_csc_mat sym_permute_lower_tri_matrix(sym_csc_mat m, i32* iperm, i32* nz_perm, sym_allocator* alloc) {
    assert(m.nrows == m.ncols);

    i32* rows = (i32*) alloc->malloc(m.nnz * sizeof(i32), alloc->ctx);
    i32* cols = (i32*) alloc->malloc(m.nnz * sizeof(i32), alloc->ctx);
    i32 col = 0;
    i32 x = iperm[0];
    for (i32 i = 0; i < m.nnz; ++i) {
        while (m.col_starts[col + 1] <= i) {
            ++col;
            x = iperm[col];
        }
        i32 row = m.row_indices[i];
        i32 y = iperm[row];

        // put the higher key first to maintain lower triangularity
        if (x > y) {
            rows[i] = x;
            cols[i] = y;
        } else {
            rows[i] = y;
            cols[i] = x;
        }
    }

    sym_csc_mat m2 = sym_csc_from_deduped_pairs(rows, cols, m.nnz, m.nrows, m.ncols, nz_perm, alloc);

    alloc->free(rows, m.nnz * sizeof(i32), alloc->ctx);
    alloc->free(cols, m.nnz * sizeof(i32), alloc->ctx);

    return m2;
}

sym_linearizer sym_linearizer_new(
    sym_csc_mat Hl_block,
    i32* Hl_block_nz_indices, i32 nblocks,
    i32* key_sizes, i32 nkeys,
    i32* key_perm,
    sym_linearization* lin,
    sym_allocator* alloc
) {
    i32* Hl_block_nz_perm = (i32*) alloc->malloc(Hl_block.nnz * sizeof(i32), alloc->ctx);
    i32* Hl_block_nz_iperm = (i32*) alloc->malloc(Hl_block.nnz * sizeof(i32), alloc->ctx);

    // TODO: take this as an argument and pipe through the one metis already computes
    i32* key_iperm = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
    for (i32 i = 0; i < nkeys; ++i) {
        key_iperm[key_perm[i]] = i;
    }

    sym_csc_mat Hl_block_permuted = sym_permute_lower_tri_matrix(Hl_block, key_iperm, Hl_block_nz_perm, alloc);

    for (i32 i = 0; i < Hl_block.nnz; ++i) {
        Hl_block_nz_iperm[Hl_block_nz_perm[i]] = i;
    }

    alloc->free(Hl_block_nz_perm, Hl_block.nnz * sizeof(i32), alloc->ctx);

    i32* key_size_scan = (i32*) alloc->malloc((nkeys + 1) * sizeof(i32), alloc->ctx);
    key_size_scan[0] = 0;
    for (i32 i = 0; i < nkeys; ++i) {
        key_size_scan[i + 1] = key_size_scan[i] + key_sizes[key_perm[i]];
    }

    i32 Hl_size = key_size_scan[nkeys];

    i32 Hl_nnz = 0;
    // TODO: rename this to reflect that it's by column?
    //   and really it's the max nnz for a triangular block
    i32* Hl_nnz_by_col_key = (i32*) alloc->malloc(nkeys * sizeof(i32), alloc->ctx);
    for (i32 i = 0; i < nkeys; ++i) {
        Hl_nnz_by_col_key[i] = 0; 
    }

    i32* Hl_key_starts_by_block_nz = (i32*) alloc->malloc(Hl_block_permuted.nnz * sizeof(i32), alloc->ctx);
    {
        i32 col = 0;
        i32 key_start = 0;
        for (i32 i = 0; i < Hl_block_permuted.nnz; ++i) {
            i32 row = Hl_block_permuted.row_indices[i];
            while (Hl_block_permuted.col_starts[col + 1] <= i) {
                ++col;
                key_start = 0;
            }
            
            i32 key_size = key_size_scan[row + 1] - key_size_scan[row];
            Hl_nnz_by_col_key[col] += key_size;

            if (row == col) {
                Hl_nnz += key_size * (key_size + 1) / 2;
            } else {
                Hl_nnz += key_size * (key_size_scan[col + 1] - key_size_scan[col]);
            }
            
            Hl_key_starts_by_block_nz[i] = key_start;
            key_start += key_size;
        }
    }

    i32* Hl_col_starts = (i32*) alloc->malloc((Hl_size + 1) * sizeof(i32), alloc->ctx);
    {
        i32 col = 0;
        Hl_col_starts[0] = 0;
        for (i32 key = 0; key < nkeys; ++key) {
            i32 key_size = key_size_scan[key + 1] - key_size_scan[key];
            for (i32 i = 0; i < key_size; ++i) {
                Hl_col_starts[col + 1] = Hl_col_starts[col] + (Hl_nnz_by_col_key[key] - i);
                ++col;
            }
        }
        assert(col == Hl_size);
        assert(Hl_col_starts[Hl_size] == Hl_nnz);
    }

    alloc->free(Hl_nnz_by_col_key, nkeys * sizeof(i32), alloc->ctx);

    i32* Hl_row_indices = (i32*) alloc->malloc(Hl_nnz * sizeof(i32), alloc->ctx);
    {
        i32 nz_index = 0;
        for (i32 col_key = 0; col_key < nkeys; ++col_key) {
            i32 col_key_size = key_size_scan[col_key + 1] - key_size_scan[col_key];
            i32 col_block_nnz = Hl_block_permuted.col_starts[col_key + 1] - Hl_block_permuted.col_starts[col_key];
            for (i32 i = 0; i < col_key_size; ++i) {
                for (i32 j = 0; j < col_block_nnz; ++j) {
                    i32 row_key = Hl_block_permuted.row_indices[Hl_block_permuted.col_starts[col_key] + j];
                    i32 row_start = key_size_scan[row_key];
                    i32 row_key_size = key_size_scan[row_key + 1] - row_start;
                    for (i32 k = (j == 0 ? i : 0); k < row_key_size; ++k) {
                        Hl_row_indices[nz_index++] = row_start + k;
                    }
                }
            }
        }
    }

    f64* Hl_data = (f64*) alloc->malloc(Hl_nnz * sizeof(f64), alloc->ctx);

    sym_csc_mat_free(Hl_block_permuted, alloc);

    sym_csc_mat Hl = {
        .nrows = Hl_size,
        .ncols = Hl_size,
        .nnz = Hl_nnz,
        .col_starts = Hl_col_starts,
        .row_indices = Hl_row_indices,
        .data = Hl_data,
    };
    lin->Hl = Hl;

    f64* rhs_data = (f64*) alloc->malloc(Hl_size * sizeof(f64), alloc->ctx);
    sym_vec rhs = {
        .n = Hl_size,
        .data = rhs_data,
    };
    lin->rhs = rhs;

    sym_linearizer lzr = {
        .nkeys = nkeys,
        .key_iperm = key_iperm,
        .key_size_scan = key_size_scan,
        .nblocks = nblocks,
        .Hl_block_nnz = Hl_block.nnz,
        .Hl_block_nz_indices = Hl_block_nz_indices,
        .Hl_key_starts_by_block_nz = Hl_key_starts_by_block_nz,
        .Hl_block_nz_iperm = Hl_block_nz_iperm,
    };
    return lzr;
}

void sym_linearizer_free(sym_linearizer lzr, sym_allocator* alloc) {
    alloc->free(lzr.key_iperm, lzr.nkeys * sizeof(i32), alloc->ctx);
    alloc->free(lzr.key_size_scan, (lzr.nkeys + 1) * sizeof(i32), alloc->ctx);
    alloc->free(lzr.Hl_key_starts_by_block_nz, lzr.Hl_block_nnz * sizeof(i32), alloc->ctx);
    alloc->free(lzr.Hl_block_nz_iperm, lzr.Hl_block_nnz * sizeof(i32), alloc->ctx);
}

// TODO: allocate lin.Hl through sym_csc_mat_new?
void sym_linearization_free(sym_linearization lin, sym_allocator* alloc) {
    sym_csc_mat_free(lin.Hl, alloc);
    sym_vec_free(lin.rhs, alloc);
}

void sym_linearization_clear(sym_linearization lin) {
    sym_csc_mat_zero(lin.Hl);
    sym_vec_zero(lin.rhs);
}

void sym_linearizer_add_hessian_tri_block(
    sym_linearizer lzr, sym_linearization lin, 
    i32 block_index, i32 key, 
    f64* data, i32 stride, i32 data_key_offset
) {
    i32 new_key = lzr.key_iperm[key];
    i32 key_start = lzr.key_size_scan[new_key];
    i32 key_size = lzr.key_size_scan[new_key + 1] - key_start;

    i32 nz_index = lzr.Hl_block_nz_indices[block_index];
    i32 new_nz_index = lzr.Hl_block_nz_iperm[nz_index];
    i32 row_offset = lzr.Hl_key_starts_by_block_nz[new_nz_index];

    for (i32 i = 0; i < key_size; ++i) {
        i32 col_start = lin.Hl.col_starts[key_start + i];
        for (i32 j = 0; j < key_size - i; ++j) {
            lin.Hl.data[col_start + (row_offset + j)] +=
                data[(data_key_offset + i) * stride + (data_key_offset + j + i)];
        }
    }
}

// for now, the row_key/col_key order must match the order in the dense block data data
void sym_linearizer_add_hessian_rect_block(
    sym_linearizer lzr, sym_linearization lin, 
    i32 block_index, i32 row_key, i32 col_key,
    f64* data, i32 stride, i32 data_row_key_offset, i32 data_col_key_offset
) {
    assert(data_row_key_offset > data_col_key_offset);

    i32 new_row_key = lzr.key_iperm[row_key];
    i32 new_col_key = lzr.key_iperm[col_key];
    i32 row_key_start = lzr.key_size_scan[new_row_key];
    i32 col_key_start = lzr.key_size_scan[new_col_key];
    i32 row_key_size = lzr.key_size_scan[new_row_key + 1] - row_key_start;
    i32 col_key_size = lzr.key_size_scan[new_col_key + 1] - col_key_start;

    i32 nz_index = lzr.Hl_block_nz_indices[block_index];
    i32 new_nz_index = lzr.Hl_block_nz_iperm[nz_index];
    i32 row_offset = lzr.Hl_key_starts_by_block_nz[new_nz_index];

    if (new_row_key > new_col_key) {
        for (i32 i = 0; i < col_key_size; ++i) {
            i32 col_start = lin.Hl.col_starts[col_key_start + i];
            for (i32 j = 0; j < row_key_size; ++j) {
                lin.Hl.data[col_start + (row_offset + j) - i] += 
                    data[(data_col_key_offset + i) * stride + (data_row_key_offset + j)];
            }
        }
    } else {
        for (i32 i = 0; i < row_key_size; ++i) {
            i32 col_start = lin.Hl.col_starts[row_key_start + i];
            for (i32 j = 0; j < col_key_size; ++j) {
                lin.Hl.data[col_start + (row_offset + j) - i] +=
                    data[(data_col_key_offset + j) * stride + (data_row_key_offset + i)];
            }
        }
    }
}

void sym_linearizer_add_rhs_block(
    sym_linearizer lzr, sym_linearization lin, i32 key, 
    f64* data, i32 data_offset) {
    i32 new_key = lzr.key_iperm[key];
    i32 key_start = lzr.key_size_scan[new_key];
    i32 key_size = lzr.key_size_scan[new_key + 1] - key_start;
    for (i32 i = 0; i < key_size; ++i) {
        lin.rhs.data[key_start + i] += data[data_offset + i];
    }
}
