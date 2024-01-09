#pragma once

#include "types.h"
#include "alloc.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    f64* data;
    i32 n;
} sym_vec;

sym_vec sym_vec_new(i32 n, sym_allocator* alloc);
void sym_vec_free(sym_vec v, sym_allocator* alloc);

void sym_vec_zero(sym_vec v);

typedef struct {
    i32* col_starts; // size ncols + 1
    i32* row_indices; // size nnz
    f64* data; // size nnz, nullable
    i32 nrows, ncols, nnz;
} sym_csc_mat;

sym_csc_mat sym_csc_mat_new(i32 nrows, i32 ncols, i32 nnz, sym_allocator* alloc);
void sym_csc_mat_free(sym_csc_mat m, sym_allocator* alloc);

void sym_csc_mat_zero(sym_csc_mat m);

void sym_sort_pairs(
    i32* rows, i32* cols, i32 n, 
    i32 nrows, i32 ncols, 
    i32* perm, 
    sym_allocator* alloc);

// nz_indices optional/nullable
sym_csc_mat sym_csc_from_pairs(
    i32* rows, i32* cols, i32 n, 
    i32 nrows, i32 ncols, 
    i32* nz_indices, 
    sym_allocator* alloc);

sym_csc_mat sym_csc_from_deduped_pairs(
    i32* rows, i32* cols, i32 n, 
    i32 nrows, i32 ncols, 
    i32* perm,
    sym_allocator* alloc
);

void sym_print_csc_mat(sym_csc_mat m, sym_allocator* alloc);


#ifdef __cplusplus
}
#endif

