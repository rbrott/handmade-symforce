#include "optimizer.h"

#include <string.h>
#include <math.h>
#include <float.h>

#include "alloc.h"
#include "linearizer.h"
#include "mat.h"
#include "sym_assert.h"


#define SWAP_PTR(Type, a, b) \
    do {                     \
        Type* tmp = a;       \
        a = b;               \
        b = tmp;             \
    } while (0)


typedef struct {
    f64 initial_lambda;
    f64 lambda_up_factor;
    f64 lambda_down_factor;
    f64 lambda_lower_bound;
    f64 lambda_upper_bound;
    f64 early_exit_min_reduction;
    f64 epsilon;
} sym_optimizer_params;

struct sym_optimizer {
    // sym_linearization* could be void* if not for the need to damp it
    // it could still be upper/lower agnostic at least
    f64 (*linearize)(sym_vec state, sym_linearization lin, void* userdata);
    void (*solve)(sym_linearization lin, sym_vec delta, void* userdata);
    void (*retract)(sym_vec state, sym_vec delta, sym_vec new_state, void* userdata);

    sym_optimizer_params params;

    f64 lambda;
    i32 iteration;

    f64 error;
    sym_vec state;
    sym_vec state_temp;
    f64* lin_rhs_temp;
    f64* lin_Hl_temp;

    sym_vec delta;

    sym_linearization lin;

    void* userdata;
};

sym_optimizer* sym_optimizer_new(
    f64 linearize(sym_vec state, sym_linearization lin, void* userdata),
    void solve(sym_linearization lin, sym_vec delta, void* userdata),
    void retract(sym_vec state, sym_vec delta, sym_vec new_state, void* userdata),
    sym_vec initial_state, sym_linearization lin,
    sym_allocator* alloc, void* userdata) {
    sym_optimizer* opt = (sym_optimizer*) alloc->malloc(sizeof(sym_optimizer), alloc->ctx);

    opt->linearize = linearize;
    opt->solve = solve;
    opt->retract = retract;

    opt->params = (sym_optimizer_params) {0};
    opt->params.initial_lambda = 1.0;
    opt->params.lambda_up_factor = 4.0;
    opt->params.lambda_down_factor = 1 / 4.0;
    opt->params.lambda_lower_bound = 0.0;
    opt->params.lambda_upper_bound = 1e10;
    opt->params.early_exit_min_reduction = 1e-6;
    opt->params.epsilon = 10.0 * DBL_EPSILON;

    opt->lambda = opt->params.initial_lambda;
    opt->iteration = 0;

    SYM_ASSERT(lin.Hl.data == NULL);
    SYM_ASSERT(lin.rhs.data == NULL);
    opt->lin = lin;
    opt->lin.Hl.data = (f64*) alloc->malloc(lin.Hl.nnz * sizeof(f64), alloc->ctx);
    opt->lin.rhs.data = (f64*) alloc->malloc(lin.rhs.n * sizeof(f64), alloc->ctx);

    opt->error = linearize(initial_state, opt->lin, userdata);
    opt->state = sym_vec_new(initial_state.n, alloc);
    memcpy(opt->state.data, initial_state.data, initial_state.n * sizeof(f64));
    opt->state_temp = sym_vec_new(initial_state.n, alloc);

    opt->lin_rhs_temp = (f64*) alloc->malloc(lin.rhs.n * sizeof(f64), alloc->ctx);
    opt->lin_Hl_temp = (f64*) alloc->malloc(lin.Hl.nnz * sizeof(f64), alloc->ctx);

    opt->delta = sym_vec_new(lin.rhs.n, alloc);

    opt->userdata = userdata;

    // this needs to compute an initial error.

    return opt;
}

sym_optimizer_status sym_optimizer_step(sym_optimizer* opt) {
    {
        // damp hessian
        sym_linearization lin_damped = opt->lin;
        lin_damped.Hl.data = opt->lin_Hl_temp;
        i32 col = 0;
        for (i32 i = 0; i < opt->lin.Hl.nnz; ++i) {
            while (opt->lin.Hl.col_starts[col + 1] <= i) {
                ++col;
            }
            if (opt->lin.Hl.row_indices[i] == col) {
                lin_damped.Hl.data[i] = opt->lin.Hl.data[i] + opt->lambda;
            } else {
                lin_damped.Hl.data[i] = opt->lin.Hl.data[i];
            }
        }

        opt->solve(lin_damped, opt->delta, opt->userdata);
    }

    opt->retract(opt->state, opt->delta, opt->state_temp, opt->userdata);

    f64 error;
    {
        sym_linearization lin_temp = opt->lin;
        lin_temp.Hl.data = opt->lin_Hl_temp;
        lin_temp.rhs.data = opt->lin_rhs_temp;
        error = opt->linearize(opt->state_temp, lin_temp, opt->userdata);
    }

    f64 relative_reduction = (opt->error - error) / (opt->error + opt->params.epsilon);

    printf("Optimizer [iter %4d] lambda: %e, error prev/new: %e/%e, rel reduction: %+e, (%a)\n",
      opt->iteration, opt->lambda, opt->error, error, relative_reduction, error);

    if (relative_reduction > -opt->params.early_exit_min_reduction / 10 &&
        relative_reduction < opt->params.early_exit_min_reduction) {
        return CONVERGED;
    }

    bool accept_update = relative_reduction > 0;

    if (!accept_update && opt->lambda >= opt->params.lambda_upper_bound) {
        return FAILED;
    }

    if (accept_update) {
        opt->lambda *= opt->params.lambda_down_factor;
        opt->error = error;

        // swap in the new buffers
        SWAP_PTR(f64, opt->lin.Hl.data, opt->lin_Hl_temp);
        SWAP_PTR(f64, opt->lin.rhs.data, opt->lin_rhs_temp);
        SWAP_PTR(f64, opt->state.data, opt->state_temp.data);
    } else {
        opt->lambda *= opt->params.lambda_up_factor;
    }

    opt->lambda = fmax(fmin(opt->lambda, opt->params.lambda_upper_bound), opt->params.lambda_lower_bound);

    ++opt->iteration;
    return IN_PROGRESS;
}

void sym_optimizer_free(sym_optimizer* opt, sym_allocator* alloc) {
    alloc->free(opt->lin_Hl_temp, opt->lin.Hl.nnz * sizeof(f64), alloc->ctx);
    alloc->free(opt->lin_rhs_temp, opt->lin.rhs.n * sizeof(f64), alloc->ctx);

    alloc->free(opt->lin.Hl.data, opt->lin.Hl.nnz * sizeof(f64), alloc->ctx);
    alloc->free(opt->lin.rhs.data, opt->lin.rhs.n * sizeof(f64), alloc->ctx);

    sym_vec_free(opt->state_temp, alloc);
    sym_vec_free(opt->state, alloc);

    sym_vec_free(opt->delta, alloc);

    alloc->free(opt, sizeof(sym_optimizer), alloc->ctx);
}
