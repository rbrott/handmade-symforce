#pragma once

#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

// see https://nullprogram.com/blog/2023/12/17/
typedef struct {
    u8 *beg;
    u8 *end;
    size nalloc, max_nalloc;
} sym_arena;

void *sym_arena_malloc(size n, void *ctx);

void sym_arena_free(void *ptr, size n, void *ctx);

#ifdef __cplusplus
}
#endif
