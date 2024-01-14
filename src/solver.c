/* ----------------------------------------------------------------------------
 * SymForce - Copyright 2022, Skydio, Inc.
 * This source code is under the MPL2 license found in the LICENSE file.
 * ---------------------------------------------------------------------------- */

// pretty much direct translation of the symforce version

#include "solver.h"

#include <assert.h>

sym_chol_solver sym_new_chol_solver(sym_csc_mat m, sym_chol_factorization* fac, sym_allocator* alloc) {
    assert(m.nrows == m.ncols);
    
    i32* visited = alloc->malloc(m.nrows * sizeof(i32), alloc->ctx);
    i32* parent = alloc->malloc(m.nrows * sizeof(i32), alloc->ctx);
    i32* nnz_by_col = alloc->malloc(m.nrows * sizeof(i32), alloc->ctx);
    for (i32 i = 0; i < m.nrows; ++i) {
        visited[i] = -1;
        parent[i] = -1;
        nnz_by_col[i] = 0;
    }

    // See Chapter 6:
    // https://www.tau.ac.il/~stoledo/Support/chapter-direct.pdf

    // Iterate through each dim k and for L(k, :), touch all nodes reachable in elimination
    // tree from nonzero entries in A(0:k-1, k)
    for (i32 k = 0; k < m.nrows; ++k) {
        // Mark k as visited
        visited[k] = k;

        // Iterate through each nonzero entry in A(0:k-1, k)
        for (i32 j = m.col_starts[k]; j < m.col_starts[k + 1]; ++j) {
            // Get row index
            i32 i = m.row_indices[j];

            // Skip if not in the upper triangle
            if (i >= k) {
                continue;
            }

            // Follow path from i to root, stop when hit previously visited node
            while (visited[i] != k) {
                // Set parent
                if (parent[i] == -1) {
                    parent[i] = k;
                }

                // L(k, i) is nonzero
                ++nnz_by_col[i];

                // Mark i as visited
                visited[i] = k;

                // Follow to parent
                i = parent[i];
            }
        }
    }

    // Allocate memory for cholesky factorization using nonzero counts
    i32* L_col_starts = alloc->malloc((m.nrows + 1) * sizeof(i32), alloc->ctx);
    L_col_starts[0] = 0;
    for (i32 i = 0; i < m.nrows; ++i) {
        L_col_starts[i + 1] = L_col_starts[i] + nnz_by_col[i];
    }
    i32 L_nnz = L_col_starts[m.nrows];
    i32* L_row_indices = alloc->malloc(L_nnz * sizeof(i32), alloc->ctx);
    f64* L_data = alloc->malloc(L_nnz * sizeof(f64), alloc->ctx);
    sym_csc_mat L = {
        .nrows = m.nrows,
        .ncols = m.nrows,
        .nnz = L_nnz,
        .col_starts = L_col_starts,
        .row_indices = L_row_indices,
        .data = L_data,
    };
    fac->L = L;

    f64* D_data = alloc->malloc(m.nrows * sizeof(f64), alloc->ctx);
    sym_vec D = {
        .n = m.nrows,
        .data = D_data,
    };
    fac->D = D;

    // Allocate other memory used for subsequent factorization and solve calls
    i32* L_k_pattern = alloc->malloc(m.nrows * sizeof(i32), alloc->ctx);
    f64* D_agg = alloc->malloc(m.nrows * sizeof(f64), alloc->ctx);

    sym_chol_solver s = {
        .visited = visited,
        .parent = parent,
        .nnz_by_col = nnz_by_col,
        .L_k_pattern = L_k_pattern,
        .D_agg = D_agg,
        .dim = m.nrows,
    };
    return s;
}

// See "Modified Cholesky Factorization", page 145:
// http://www.bioinfo.org.cn/~wangchao/maa/Numerical_Optimization.pdf

// TODO: row indices is computed each time -- maybe that's fine?
void sym_chol_solver_factor(sym_chol_solver solver, sym_csc_mat m, sym_chol_factorization fac) {
    assert(fac.L.nrows == fac.L.ncols);
    assert(solver.dim == fac.L.nrows);

    // Initialize helpers
    for (i32 i = 0; i < solver.dim; ++i) {
        solver.nnz_by_col[i] = 0;
        solver.D_agg[i] = 0.0;
    }

    // For each row of L, compute nonzero pattern in topo order
    for (i32 k = 0; k < solver.dim; ++k) {
        // Mark k as visited
        solver.visited[k] = k;

        // Reverse counter
        i32 top_inx = solver.dim;

        for (i32 j = m.col_starts[k]; j < m.col_starts[k + 1]; ++j) {
            // Get row index
            i32 i = m.row_indices[j];

            // Skip if not in the upper triangle
            if (i > k) {
                continue;
            }

            // Sum A(i, k) into D_agg
            solver.D_agg[i] += m.data[j];

            // Follow path from i to root, stop when hit previously visited node
            i32 depth = 0;
            while (solver.visited[i] != k) {
                // L(k,i) is nonzero
                solver.L_k_pattern[depth] = i;

                // Mark i as visited
                solver.visited[i] = k;

                // Follow to parent
                i = solver.parent[i];

                // Increment depth
                ++depth;
            }

            // Update pattern
            while (depth > 0) {
                --depth;
                --top_inx;
                solver.L_k_pattern[top_inx] = solver.L_k_pattern[depth];
            }
        }

        // Get D(k, k) and clear D_agg(k)
        f64 D_k = solver.D_agg[k];
        solver.D_agg[k] = 0.0;

        // NOTE: This is a double loop in a loop and is ~O(N^3 / 6)
        for (; top_inx < solver.dim; ++top_inx) {
            // L_k_pattern_[top_inx:] is the pattern of L(:, k)
            i32 i = solver.L_k_pattern[top_inx];

            // Compute the nonzero L(k, i)
            f64 D_agg_i = solver.D_agg[i];
            f64 L_ki = D_agg_i / fac.D.data[i];

            // Get the range for i
            i32 ptr_start = fac.L.col_starts[i];
            i32 ptr_end = ptr_start + solver.nnz_by_col[i];

            // Update D_agg
            solver.D_agg[i] = 0.0;
            for (i32 ptr = ptr_start; ptr < ptr_end; ++ptr) {
                solver.D_agg[fac.L.row_indices[ptr]] -= fac.L.data[ptr] * D_agg_i;
            }

            // Save L(k, i)
            fac.L.row_indices[ptr_end] = k;
            fac.L.data[ptr_end] = L_ki;

            // Update D(k)
            D_k -= L_ki * D_agg_i;

            // Increment nonzeros in column i
            ++solver.nnz_by_col[i];
        }

        // Save D(k)
        fac.D.data[k] = D_k;
    }
}

// TODO: allocation here is pretty sad
//   not sure if making transpose not alloc is feasible
//   I can't tell what Eigen does
void sym_chol_solver_solve_in_place(sym_chol_factorization fac, sym_vec x, sym_allocator* alloc) {
    assert(fac.L.nrows == fac.L.ncols);
    assert(x.n == fac.L.nrows);

    // L \ x in place
    for (i32 k = 0; k < fac.L.nrows; ++k) {
        i32 j = fac.L.col_starts[k];
        // TODO: apparently the diagonal is all empty?
        // assert(fac.L.row_indices[j] == k);
        f64 y = x.data[k];
        for (; j < fac.L.col_starts[k + 1]; ++j) {
            i32 i = fac.L.row_indices[j];
            x.data[i] -= fac.L.data[j] * y;
        }
    }

    // D \ x in place
    for (i32 i = 0; i < fac.D.n; ++i) {
        x.data[i] /= fac.D.data[i];
    }

    // L^T \ x in place
    i32* perm = alloc->malloc(fac.L.nnz * sizeof(i32), alloc->ctx);
    sym_csc_mat Lt = sym_transpose_csc(fac.L, perm, alloc);
    alloc->free(perm, fac.L.nnz * sizeof(i32), alloc->ctx);
    for (i32 k = Lt.nrows - 1; k >= 0; --k) {
        // make sure there's a diagonal entry
        i32 j = Lt.col_starts[k + 1] - 1;
        f64 y = x.data[k];
        for (; j >= Lt.col_starts[k]; --j) {
            i32 i = Lt.row_indices[j];
            x.data[i] -= Lt.data[j] * y;
        }
    }
        
    sym_csc_mat_free(Lt, alloc);
}

void sym_chol_solver_free(sym_chol_solver solver, sym_allocator* alloc) {
    alloc->free(solver.visited, solver.dim * sizeof(i32), alloc->ctx);
    alloc->free(solver.parent, solver.dim * sizeof(i32), alloc->ctx);
    alloc->free(solver.nnz_by_col, solver.dim * sizeof(i32), alloc->ctx);
    alloc->free(solver.L_k_pattern, solver.dim * sizeof(i32), alloc->ctx);
    alloc->free(solver.D_agg, solver.dim * sizeof(f64), alloc->ctx);
}

void sym_chol_factorization_free(sym_chol_factorization fac, sym_allocator* alloc) {
    sym_csc_mat_free(fac.L, alloc);
    sym_vec_free(fac.D, alloc);
}
