#include "timing.h"

#include <stdlib.h>
#include <string.h>

#include "sym_assert.h"

#if defined(_WIN32)
#include <windows.h>
#else
#include <time.h>
#endif

enum {
    SYM_TIME_INITIAL_CAPACITY = 16,
};

typedef struct {
    char *name;
    sym_time_stats stats;
} sym_time_record;

struct sym_timing {
    sym_allocator *alloc;
    sym_time_record *records;
    i32 count;
    i32 capacity;
};

static u64 sym_time_now_ns(void) {
#if defined(_WIN32)
    LARGE_INTEGER counter;
    LARGE_INTEGER frequency;
    QueryPerformanceCounter(&counter);
    QueryPerformanceFrequency(&frequency);

    const u64 nanoseconds_per_second = 1000000000;
    const u64 seconds = (u64) counter.QuadPart / (u64) frequency.QuadPart;
    const u64 remainder = (u64) counter.QuadPart % (u64) frequency.QuadPart;
    return seconds * nanoseconds_per_second
        + remainder * nanoseconds_per_second / (u64) frequency.QuadPart;
#else
    struct timespec now;
    const int result = clock_gettime(CLOCK_MONOTONIC, &now);
    SYM_ASSERT(result == 0);

    const u64 nanoseconds_per_second = 1000000000;
    return (u64) now.tv_sec * nanoseconds_per_second + (u64) now.tv_nsec;
#endif
}

static sym_time_record *sym_time_find(
    const sym_timing *timing,
    const char *name
) {
    for (i32 i = 0; i < timing->count; ++i) {
        if (strcmp(timing->records[i].name, name) == 0) {
            return timing->records + i;
        }
    }

    return NULL;
}

static void sym_time_grow(sym_timing *timing) {
    const i32 capacity = timing->capacity * 2;
    const size bytes = (size) capacity * sizeof(sym_time_record);
    sym_time_record *records = timing->alloc->malloc(bytes, timing->alloc->ctx);

    memcpy(
        records,
        timing->records,
        (usize) timing->count * sizeof(sym_time_record)
    );

    const size old_bytes = (size) timing->capacity * sizeof(sym_time_record);
    timing->alloc->free(timing->records, old_bytes, timing->alloc->ctx);

    timing->records = records;
    timing->capacity = capacity;
}

static sym_time_record *sym_time_add(sym_timing *timing, const char *name) {
    if (timing->count == timing->capacity) {
        sym_time_grow(timing);
    }

    const usize name_size = strlen(name) + 1;
    char *name_copy = timing->alloc->malloc((size) name_size, timing->alloc->ctx);
    memcpy(name_copy, name, name_size);

    sym_time_record *record = timing->records + timing->count;
    *record = (sym_time_record) {
        .name = name_copy,
        .stats = {
            .min_ns = UINT64_MAX,
        },
    };
    ++timing->count;

    return record;
}

static int sym_time_compare(const void *left, const void *right) {
    const sym_time_record *const *a = left;
    const sym_time_record *const *b = right;

    if ((*a)->stats.total_ns < (*b)->stats.total_ns) {
        return 1;
    }
    if ((*a)->stats.total_ns > (*b)->stats.total_ns) {
        return -1;
    }
    return strcmp((*a)->name, (*b)->name);
}

sym_timing *sym_timing_new(sym_allocator *alloc) {
    sym_timing *timing = alloc->malloc(sizeof(sym_timing), alloc->ctx);
    const size records_size = SYM_TIME_INITIAL_CAPACITY * sizeof(sym_time_record);

    *timing = (sym_timing) {
        .alloc = alloc,
        .records = alloc->malloc(records_size, alloc->ctx),
        .capacity = SYM_TIME_INITIAL_CAPACITY,
    };

    return timing;
}

void sym_timing_free(sym_timing *timing) {
    if (timing == NULL) {
        return;
    }

    sym_allocator *alloc = timing->alloc;
    sym_timing_reset(timing);

    const size records_size = (size) timing->capacity * sizeof(sym_time_record);
    alloc->free(timing->records, records_size, alloc->ctx);
    alloc->free(timing, sizeof(sym_timing), alloc->ctx);
}

void sym_timing_reset(sym_timing *timing) {
    for (i32 i = 0; i < timing->count; ++i) {
        const size name_size = (size) strlen(timing->records[i].name) + 1;
        timing->alloc->free(timing->records[i].name, name_size, timing->alloc->ctx);
    }

    timing->count = 0;
}

void sym_timing_record(sym_timing *timing, const char *name, u64 duration_ns) {
    sym_time_record *record = sym_time_find(timing, name);
    if (record == NULL) {
        record = sym_time_add(timing, name);
    }

    sym_time_stats *stats = &record->stats;
    ++stats->count;
    stats->total_ns += duration_ns;
    stats->min_ns = duration_ns < stats->min_ns ? duration_ns : stats->min_ns;
    stats->max_ns = duration_ns > stats->max_ns ? duration_ns : stats->max_ns;
}

b32 sym_timing_get(
    const sym_timing *timing,
    const char *name,
    sym_time_stats *stats
) {
    const sym_time_record *record = sym_time_find(timing, name);
    if (record == NULL) {
        return 0;
    }

    *stats = record->stats;
    return 1;
}

void sym_timing_print(const sym_timing *timing, FILE *file) {
    if (timing->count == 0) {
        return;
    }

    const size sorted_size = (size) timing->count * sizeof(sym_time_record *);
    sym_time_record **sorted = timing->alloc->malloc(sorted_size, timing->alloc->ctx);
    for (i32 i = 0; i < timing->count; ++i) {
        sorted[i] = timing->records + i;
    }
    qsort(sorted, (usize) timing->count, sizeof(sym_time_record *), sym_time_compare);

    const f64 nanoseconds_per_second = 1000000000.0;
    fprintf(file, "\nTiming results:\n");
    fprintf(file, "%-28s %8s %12s %12s %12s %12s\n",
        "Name", "Count", "Total (s)", "Mean (s)", "Min (s)", "Max (s)");

    for (i32 i = 0; i < timing->count; ++i) {
        const sym_time_record *record = sorted[i];
        const sym_time_stats stats = record->stats;
        const f64 total = (f64) stats.total_ns / nanoseconds_per_second;
        const f64 mean = total / (f64) stats.count;
        const f64 min = (f64) stats.min_ns / nanoseconds_per_second;
        const f64 max = (f64) stats.max_ns / nanoseconds_per_second;

        fprintf(file, "%-28s %8llu %12.6f %12.6f %12.6f %12.6f\n",
            record->name,
            (unsigned long long) stats.count,
            total,
            mean,
            min,
            max);
    }

    timing->alloc->free(sorted, sorted_size, timing->alloc->ctx);
}

sym_time_scope sym_tic(sym_timing *timing, const char *name) {
    return (sym_time_scope) {
        .timing = timing,
        .name = name,
        .start_ns = sym_time_now_ns(),
        .active = 1,
    };
}

void sym_toc(sym_time_scope *scope) {
    if (!scope->active) {
        return;
    }

    const u64 duration_ns = sym_time_now_ns() - scope->start_ns;
    sym_timing_record(scope->timing, scope->name, duration_ns);
    scope->active = 0;
}
