#pragma once

#include "alloc.h"
#include "linearizer.h"
#include "mat.h"

typedef enum {
    INVALID = 0,
    CONVERGED = 1,
    IN_PROGRESS = 2,
    FAILED = 3,
} sym_optimizer_status;

typedef struct sym_optimizer sym_optimizer;

sym_optimizer* sym_optimizer_new(
    // linearizer returns the error
    f64 linearize(sym_vec state, sym_linearization lin, void* userdata),
    void solve(sym_linearization lin, sym_vec delta, void* userdata),
    void retract(sym_vec state, sym_vec delta, sym_vec new_state, void* userdata),
    sym_vec initial_state, sym_linearization lin,
    sym_allocator* alloc, void* userdata);

void sym_optimizer_free(sym_optimizer* opt, sym_allocator* alloc);

sym_optimizer_status sym_optimizer_step(sym_optimizer* opt);
