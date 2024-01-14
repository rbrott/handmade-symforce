/* ----------------------------------------------------------------------------
 * SymForce - Copyright 2022, Skydio, Inc.
 * This source code is under the MPL2 license found in the LICENSE file.
 * ---------------------------------------------------------------------------- */

// pretty much direct translation of the symforce version

#include "types.h"
#include "mat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    i32* visited; // size dim
    i32* parent; // size dim
    i32* nnz_by_col; // size dim
    i32* L_k_pattern; // size dim
    f64* D_agg; // size dim
    i32* Lt_perm; // size L_nnz
    i32 dim;
    i32 L_nnz;
} sym_chol_solver;

typedef struct {
    sym_csc_mat L;
    sym_vec D;
    sym_csc_mat Lt;
} sym_chol_factorization;

sym_chol_solver sym_new_chol_solver(sym_csc_mat m, sym_chol_factorization* fac, sym_allocator* alloc);

void sym_chol_solver_factor(sym_chol_solver solver, sym_csc_mat m, sym_chol_factorization fac);

void sym_chol_solver_solve_in_place(sym_chol_factorization fac, sym_vec x);

void sym_chol_solver_free(sym_chol_solver solver, sym_allocator* alloc);

void sym_chol_factorization_free(sym_chol_factorization fac, sym_allocator* alloc);

#ifdef __cplusplus
}
#endif
