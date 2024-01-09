#pragma once

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

// taken from https://nullprogram.com/blog/2023/12/17/
// implementations need not be thread safe
typedef struct {
    // TODO: guarantee that allocations are zeroed?
    void *(*malloc)(size, void *ctx);
    void  (*free)(void *, size, void *ctx);
    void   *ctx;
} sym_allocator;

#ifdef __cplusplus
}
#endif
