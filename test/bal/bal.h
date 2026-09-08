#pragma once

#include "types.h"
#include "linearizer.h"

typedef struct {
    i32* camera_indices;
    i32* point_indices;
    f64* pixels;
    f64* values;
    i32 num_cameras;
    i32 num_points;
    i32 num_observations;
} bal_problem;

typedef struct {
    bal_problem problem;

    sym_linearizer lzr;

    void* userdata;
} bal_optimizer_state;

bal_problem bal_read_new(char* filepath, sym_allocator* alloc);
void bal_free(bal_problem p, sym_allocator* alloc);

sym_linearizer bal_linearizer_new(bal_problem p, sym_linearization* lin,
    i32** Hl_block_nz_indices, sym_allocator* alloc);

// userdata must be bal_optimizer_state*
f64 bal_linearize(sym_vec state, sym_linearization lin, void* userdata);
void bal_retract(sym_vec state, sym_vec delta, sym_vec new_state, void* userdata);
