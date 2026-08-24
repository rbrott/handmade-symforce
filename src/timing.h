#pragma once

#include <stdio.h>

#include "alloc.h"
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct sym_timing sym_timing;

typedef struct {
    u64 count;
    u64 total_ns;
    u64 min_ns;
    u64 max_ns;
} sym_time_stats;

typedef struct {
    sym_timing *timing;
    const char *name;
    u64 start_ns;
    b32 active;
} sym_time_scope;

sym_timing *sym_timing_new(sym_allocator *alloc);

void sym_timing_free(sym_timing *timing);

void sym_timing_reset(sym_timing *timing);

void sym_timing_record(sym_timing *timing, const char *name, u64 duration_ns);

b32 sym_timing_get(
    const sym_timing *timing,
    const char *name,
    sym_time_stats *stats
);

void sym_timing_print(const sym_timing *timing, FILE *file);

sym_time_scope sym_tic(sym_timing *timing, const char *name);

void sym_toc(sym_time_scope *scope);

#ifdef __cplusplus
}
#endif

#if defined(SYM_ENABLE_TIMING)
#if !defined(__clang__) && !defined(__GNUC__)
#error SYM_TIME_SCOPE requires cleanup attribute support
#endif

#define SYM_TIME_JOIN1(a, b) a##b
#define SYM_TIME_JOIN(a, b) SYM_TIME_JOIN1(a, b)
#define SYM_TIME_SCOPE(timing, name)                                      \
    sym_time_scope SYM_TIME_JOIN(sym_time_scope_, __LINE__)               \
        __attribute__((cleanup(sym_toc))) = sym_tic((timing), (name))
#else
#define SYM_TIME_SCOPE(timing, name)
#endif
