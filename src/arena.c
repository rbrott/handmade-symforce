#include "arena.h"

#include <assert.h>
#include <stdio.h>

void *sym_arena_malloc(size n, void *ctx) {
    sym_arena *a = ctx;
    size available = a->end - a->beg;
    // TODO: custom alignment?
    size alignment = -n & 15;
    assert(n <= available - alignment);
    a->nalloc += n;
    return a->end -= n + alignment;
}

void sym_arena_free(void *ptr, size n, void *ctx) {
    sym_arena *a = ctx;

    a->max_nalloc = a->nalloc > a->max_nalloc ? a->nalloc : a->max_nalloc;
    a->nalloc -= n;
}
