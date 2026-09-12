/* ----------------------------------------------------------------------------
 * SymForce - Copyright 2022, Skydio, Inc.
 * This source code is under the Apache 2.0 license found in the LICENSE file.
 * ---------------------------------------------------------------------------- */

#include <float.h>
#include <stdbool.h>
#include <string.h>

#include <cholmod.h>

#include "optimizer.h"
#include "sym_assert.h"
#include "alloc.h"
#include "arena.h"
#include "linearizer.h"
#include "solver.h"
#include "types.h"
#include "bal.h"

typedef struct {
    cholmod_common chol_common;
    cholmod_sparse H;
    cholmod_factor* H_fac;
    cholmod_dense b;
    cholmod_dense* solve_y;
    cholmod_dense* solve_e;
} bal_solver_state;

void bal_solve(
    sym_linearization lin,
    sym_vec delta,
    void* userdata
) {
    bal_optimizer_state* opt_state = (bal_optimizer_state*) userdata;
    bal_solver_state* solver_state = (bal_solver_state*) opt_state->userdata;

    solver_state->H.xtype = CHOLMOD_REAL;
    solver_state->H.x = lin.Hl.data;

    int success = cholmod_factorize(&solver_state->H, solver_state->H_fac, &solver_state->chol_common);
    SYM_ASSERT(success); // TODO: idk what the value of this is supposed to be

    for (i32 i = 0; i < lin.Hl.nrows; ++i) {
        ((f64*) solver_state->b.x)[i] = -lin.rhs.data[i];
    }

    cholmod_dense x = {
        .nrow = delta.n,
        .ncol = 1,
        .nzmax = delta.n,
        .d = delta.n,
        .x = delta.data,
        .z = NULL,
        .xtype = CHOLMOD_REAL,
        .dtype = CHOLMOD_DOUBLE
    };
    cholmod_dense* xref = &x;

    success = cholmod_solve2(
        CHOLMOD_A,
        solver_state->H_fac,
        &solver_state->b,
        NULL,
        &xref,
        NULL,
        &solver_state->solve_y,
        &solver_state->solve_e,
        &solver_state->chol_common
    );
    SYM_ASSERT(success);
}

void print_error(int status, const char *file, int line,
      const char *message) {
    printf("status: %d, file: %s, line: %d, message: %s\n", status, file, line, message);
}

int main(int argc, char** argv) {
    SYM_ASSERT(argc == 2 || argc == 3);

    bool use_supernodal = argc == 3;
    if (use_supernodal) {
        SYM_ASSERT(strcmp(argv[2], "--supernodal") == 0);
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

    bal_solver_state bal_solver_state = {0};

    cholmod_start(&bal_solver_state.chol_common);
    bal_solver_state.chol_common.nmethods = 1;
    bal_solver_state.chol_common.method[0].ordering = CHOLMOD_NATURAL;
    bal_solver_state.chol_common.postorder = true;
    bal_solver_state.chol_common.supernodal = use_supernodal ? CHOLMOD_SUPERNODAL : CHOLMOD_SIMPLICIAL;

    bal_solver_state.chol_common.error_handler = print_error;

    bal_solver_state.H = (cholmod_sparse) {
        .nrow = lin.Hl.nrows,
        .ncol = lin.Hl.ncols,
        .nzmax = lin.Hl.nnz,
        .p = lin.Hl.col_starts,
        .i = lin.Hl.row_indices,
        .nz = NULL,
        .x = NULL,
        .z = NULL,
        .stype = -1, // symmetric, lower part stored
        .itype = CHOLMOD_INT,
        .xtype = CHOLMOD_PATTERN,
        .dtype = CHOLMOD_DOUBLE,
        .sorted = 1,
        .packed = 1,
    };

    // H, H_fac symbolic for now -- does this allocate memory for the data members?
    bal_solver_state.H_fac = cholmod_analyze(&bal_solver_state.H, &bal_solver_state.chol_common);

    bal_solver_state.b = (cholmod_dense) {
        .nrow = lin.Hl.nrows,
        .ncol = 1,
        .nzmax = lin.Hl.nrows,
        .d = lin.Hl.nrows,
        .x = alloc->malloc(lin.Hl.nrows * sizeof(f64), alloc->ctx),
        .z = NULL,
        .xtype = CHOLMOD_REAL,
        .dtype = CHOLMOD_DOUBLE
    };
    bal_solver_state.solve_y = NULL;
    bal_solver_state.solve_e = NULL;

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

    cholmod_free_dense(&bal_solver_state.solve_y, &bal_solver_state.chol_common);
    cholmod_free_dense(&bal_solver_state.solve_e, &bal_solver_state.chol_common);
    cholmod_free_factor(&bal_solver_state.H_fac, &bal_solver_state.chol_common);
    cholmod_finish(&bal_solver_state.chol_common);

    sym_linearizer_free(lzr, alloc);
    sym_linearization_free(lin, alloc);

    alloc->free(Hl_block_nz_indices, lzr.nblocks * sizeof(i32), alloc->ctx);

    bal_free(p, alloc);

    printf("nalloc = %ld\n", arena.nalloc);
    printf("max_nalloc = %ld\n", arena.max_nalloc);
    SYM_ASSERT(arena.nalloc == 0);

    free(buf);
}
