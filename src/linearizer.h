#pragma once

#include "types.h"
#include "alloc.h"
#include "mat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    sym_csc_mat Hl;
    sym_vec rhs;
} sym_linearization;

typedef struct {
    i32* key_iperm; // size nkeys
    i32* key_size_scan; // size nkeys + 1, post-permutation
    i32* Hl_block_nz_indices; // size nblocks, unowned
    i32* Hl_key_starts_by_block_nz; // size Hl_block_nnz
    i32* Hl_block_nz_iperm; // size Hl_block_nnz
    i32 nkeys;
    i32 nblocks;
    i32 Hl_block_nnz;
} sym_linearizer;

void sym_get_metis_tri_perm(sym_csc_mat m, i32* weights, i32* perm, sym_allocator* alloc);

// Hl_block_nz_indices must outlive the linearizer
sym_linearizer sym_linearizer_new(
    sym_csc_mat Hl_block,
    i32* Hl_block_nz_indices, i32 nblocks,
    i32* key_sizes, i32 nkeys,
    i32* key_perm, 
    sym_linearization* lin,
    sym_allocator* alloc);

void sym_linearization_clear(sym_linearization lin);

// data is a dense, column-major matrix with column size stride
void sym_linearizer_add_hessian_tri_block(
    sym_linearizer lzr, sym_linearization lin, 
    i32 block_index, i32 key, 
    f64* data, i32 stride, i32 data_key_offset);
void sym_linearizer_add_hessian_rect_block(
    sym_linearizer lzr, sym_linearization lin, 
    i32 block_index, i32 row_key, i32 col_key,
    f64* data, i32 stride, i32 data_row_key_offset, i32 data_col_key_offset);

void sym_linearizer_add_rhs_block(
    sym_linearizer lzr, sym_linearization lin, i32 key, 
    f64* data, i32 data_offset);

void sym_linearizer_free(sym_linearizer lzr, sym_allocator* alloc);
void sym_linearization_free(sym_linearization lin, sym_allocator* alloc);

#ifdef __cplusplus
}
#endif
