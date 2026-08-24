#include "timing.h"

#include <assert.h>
#include <stdlib.h>

#include "arena.h"

static void early_return(sym_timing *timing) {
    SYM_TIME_SCOPE(timing, "scope");
    return;
}

int main(void) {
    const size arena_size = 4096;
    u8 *buffer = malloc(arena_size);
    sym_arena arena = {
        .beg = buffer,
        .end = buffer + arena_size,
    };
    sym_allocator alloc = {
        .malloc = sym_arena_malloc,
        .free = sym_arena_free,
        .ctx = &arena,
    };
    sym_timing *timing = sym_timing_new(&alloc);

    sym_timing_record(timing, "alpha", 10);
    sym_timing_record(timing, "alpha", 30);
    sym_timing_record(timing, "beta", 5);

    sym_time_stats stats;
    assert(sym_timing_get(timing, "alpha", &stats));
    assert(stats.count == 2);
    assert(stats.total_ns == 40);
    assert(stats.min_ns == 10);
    assert(stats.max_ns == 30);
    assert(!sym_timing_get(timing, "missing", &stats));

    early_return(timing);
    assert(sym_timing_get(timing, "scope", &stats));
    assert(stats.count == 1);

    sym_timing_reset(timing);
    assert(!sym_timing_get(timing, "alpha", &stats));

    sym_timing_free(timing);
    assert(arena.nalloc == 0);
    free(buffer);
}
