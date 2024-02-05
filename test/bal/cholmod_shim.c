#include "sym_assert.h"

#include <cholmod_internal.h>

int CHOLMOD(amd)
(
    // input:
    cholmod_sparse *A,  // matrix to order
    Int *fset,          // subset of 0:(A->ncol)-1
    size_t fsize,       // size of fset
    // output:
    Int *Perm,          // size A->nrow, output permutation
    cholmod_common *Common
) {
    SYM_ASSERT(0);
}

int CHOLMOD(colamd)
(
    // input:
    cholmod_sparse *A,  // matrix to order
    Int *fset,          // subset of 0:(A->ncol)-1
    size_t fsize,       // size of fset
    int postorder,      // if TRUE, follow with a coletree postorder
    // output:
    Int *Perm,          // size A->nrow, output permutation
    cholmod_common *Common
)
{
    SYM_ASSERT(0);
}
