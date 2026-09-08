/* ----------------------------------------------------------------------------
 * SymForce - Copyright 2022, Skydio, Inc.
 * This source code is under the Apache 2.0 license found in the LICENSE file.
 * ---------------------------------------------------------------------------- */

#include <float.h>
#include <stdbool.h>
#include <string.h>

#include "optimizer.h"
#include "sym_assert.h"
#include "alloc.h"
#include "arena.h"
#include "linearizer.h"
#include "solver.h"
#include "types.h"
#include "bal.h"

typedef struct {
    i32* Hlt_perm;
    sym_csc_mat Hlt;

    sym_chol_factorization fac;
    sym_chol_solver solver;
} bal_solver_state;

void bal_solve(
    sym_linearization lin,
    sym_vec delta,
    void* userdata
) {
    bal_optimizer_state* opt_state = (bal_optimizer_state*) userdata;
    bal_solver_state* solver_state = (bal_solver_state*) opt_state->userdata;

    for (i32 i = 0; i < lin.Hl.nnz; ++i) {
        solver_state->Hlt.data[i] = lin.Hl.data[solver_state->Hlt_perm[i]];
    }

    sym_chol_solver_factor(solver_state->solver, solver_state->Hlt, solver_state->fac);

    for (i32 i = 0; i < lin.Hl.nrows; ++i) {
        delta.data[i] = -lin.rhs.data[i];
    }

    sym_chol_solver_solve_in_place(solver_state->fac, delta);
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

    sym_linearization lin;
    i32* Hl_block_nz_indices;
    sym_linearizer lzr = bal_linearizer_new(p, &lin, &Hl_block_nz_indices, alloc);

    i32* Hlt_perm = (i32*) alloc->malloc(lin.Hl.nnz * sizeof(i32), alloc->ctx);
    sym_csc_mat Hlt = sym_transpose_csc(lin.Hl, Hlt_perm, alloc);
    Hlt.data = (f64*) alloc->malloc(lin.Hl.nnz * sizeof(f64), alloc->ctx);

    sym_chol_factorization fac = {};
    sym_chol_solver solver = sym_new_chol_solver(Hlt, &fac, populate_lt, alloc);

    bal_solver_state bal_solver_state = {
        .Hlt_perm = Hlt_perm,
        .Hlt = Hlt,
        .fac = fac,
        .solver = solver
    };

    bal_optimizer_state bal_opt_state = {
        .problem = p,
        .lzr = lzr,
        .userdata = &bal_solver_state
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

    alloc->free(Hl_block_nz_indices, lzr.nblocks * sizeof(i32), alloc->ctx);

    bal_free(p, alloc);

    printf("nalloc = %ld\n", arena.nalloc);
    printf("max_nalloc = %ld\n", arena.max_nalloc);
    SYM_ASSERT(arena.nalloc == 0);

    free(buf);
}
