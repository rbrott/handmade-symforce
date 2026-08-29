const std = @import("std");
const zcc = @import("compile_commands");

const eigen_dir = "third_party/eigen-5.0.1";
const metis_dir = "third_party/metis-5.1.0";
const suitesparse_dir = "third_party/suite-sparse-7.6.0";

// Run `python_exe -c <expr>` at configure time and return its stdout.
fn pythonQuery(b: *std.Build, python_exe: []const u8, expr: []const u8) []const u8 {
    const result = std.process.run(b.allocator, b.graph.io, .{
        .argv = &.{ python_exe, "-c", expr },
    }) catch @panic("Missing python");
    b.allocator.free(result.stderr);
    return result.stdout;
}

// GKlib source files (relative to metis-5.1.0/GKlib).
const gklib_files = [_][]const u8{
    "htable.c",
    "itemsets.c",
    "memory.c",
    "gkregex.c",
    "util.c",
    "omp.c",
    "io.c",
    "fkvkselect.c",
    "sort.c",
    "evaluate.c",
    "pqueue.c",
    "string.c",
    "random.c",
    "fs.c",
    "pdb.c",
    "csr.c",
    "timers.c",
    "error.c",
    "seq.c",
    "b64.c",
    "mcore.c",
    "graph.c",
    "blas.c",
    "getopt.c",
    "tokenizer.c",
    "rw.c",
};

// libmetis source files (relative to metis-5.1.0/libmetis).
const metis_files = [_][]const u8{
    "kwayrefine.c",
    "mincover.c",
    "bucketsort.c",
    "parmetis.c",
    "util.c",
    "kmetis.c",
    "meshpart.c",
    "compress.c",
    "gklib.c",
    "auxapi.c",
    "separator.c",
    "frename.c",
    "mcutil.c",
    "ometis.c",
    "wspace.c",
    "sfm.c",
    "debug.c",
    "balance.c",
    "pmetis.c",
    "mmd.c",
    "refine.c",
    "contig.c",
    "coarsen.c",
    "kwayfm.c",
    "stat.c",
    "checkgraph.c",
    "timing.c",
    "fm.c",
    "fortran.c",
    "initpart.c",
    "graph.c",
    "srefine.c",
    "mesh.c",
    "options.c",
    "minconn.c",
};

// CHOLMOD source files (relative to suite-sparse-7.6.0).
const cholmod_files = [_][]const u8{
    "SuiteSparse_config/SuiteSparse_config.c",
    "CHOLMOD/Cholesky/cholmod_factorize.c",
    "CHOLMOD/Cholesky/cholmod_postorder.c",
    "CHOLMOD/Cholesky/cholmod_rowcolcounts.c",
    "CHOLMOD/Cholesky/cholmod_l_rowfac.c",
    "CHOLMOD/Cholesky/cholmod_rcond.c",
    "CHOLMOD/Cholesky/cholmod_resymbol.c",
    "CHOLMOD/Cholesky/cholmod_etree.c",
    "CHOLMOD/Cholesky/cholmod_l_solve.c",
    "CHOLMOD/Cholesky/cholmod_solve.c",
    "CHOLMOD/Cholesky/cholmod_rowfac.c",
    "CHOLMOD/Cholesky/cholmod_l_rowcolcounts.c",
    "CHOLMOD/Cholesky/cholmod_l_postorder.c",
    "CHOLMOD/Cholesky/cholmod_l_factorize.c",
    "CHOLMOD/Cholesky/cholmod_l_etree.c",
    "CHOLMOD/Cholesky/cholmod_analyze.c",
    "CHOLMOD/Cholesky/cholmod_l_spsolve.c",
    "CHOLMOD/Cholesky/cholmod_l_rcond.c",
    "CHOLMOD/Cholesky/cholmod_l_analyze.c",
    "CHOLMOD/Cholesky/cholmod_spsolve.c",
    "CHOLMOD/Cholesky/cholmod_l_resymbol.c",
    "CHOLMOD/Utility/cholmod_l_hypot.c",
    "CHOLMOD/Utility/cholmod_l_aat.c",
    "CHOLMOD/Utility/cholmod_defaults.c",
    "CHOLMOD/Utility/cholmod_xtype.c",
    "CHOLMOD/Utility/cholmod_l_transpose_unsym.c",
    "CHOLMOD/Utility/cholmod_ensure_dense.c",
    "CHOLMOD/Utility/cholmod_l_mult_size_t.c",
    "CHOLMOD/Utility/cholmod_score_comp.c",
    "CHOLMOD/Utility/cholmod_l_sbound.c",
    "CHOLMOD/Utility/cholmod_l_dense_to_sparse.c",
    "CHOLMOD/Utility/cholmod_allocate_triplet.c",
    "CHOLMOD/Utility/cholmod_l_copy_triplet.c",
    "CHOLMOD/Utility/cholmod_set_empty.c",
    "CHOLMOD/Utility/cholmod_dbound.c",
    "CHOLMOD/Utility/cholmod_sort.c",
    "CHOLMOD/Utility/cholmod_allocate_factor.c",
    "CHOLMOD/Utility/cholmod_realloc_multiple.c",
    "CHOLMOD/Utility/cholmod_l_reallocate_column.c",
    "CHOLMOD/Utility/cholmod_l_start.c",
    "CHOLMOD/Utility/cholmod_allocate_dense.c",
    "CHOLMOD/Utility/cholmod_l_alloc_work.c",
    "CHOLMOD/Utility/cholmod_nnz.c",
    "CHOLMOD/Utility/cholmod_reallocate_triplet.c",
    "CHOLMOD/Utility/cholmod_reallocate_column.c",
    "CHOLMOD/Utility/cholmod_l_copy_dense.c",
    "CHOLMOD/Utility/cholmod_add.c",
    "CHOLMOD/Utility/cholmod_change_factor.c",
    "CHOLMOD/Utility/cholmod_reallocate_factor.c",
    "CHOLMOD/Utility/cholmod_copy.c",
    "CHOLMOD/Utility/cholmod_allocate_work.c",
    "CHOLMOD/Utility/cholmod_l_free_work.c",
    "CHOLMOD/Utility/cholmod_l_copy_dense2.c",
    "CHOLMOD/Utility/cholmod_l_reallocate_factor.c",
    "CHOLMOD/Utility/cholmod_pack_factor.c",
    "CHOLMOD/Utility/cholmod_l_ones.c",
    "CHOLMOD/Utility/cholmod_l_dense_nnz.c",
    "CHOLMOD/Utility/cholmod_l_calloc.c",
    "CHOLMOD/Utility/cholmod_mult_uint64_t.c",
    "CHOLMOD/Utility/cholmod_alloc_factor.c",
    "CHOLMOD/Utility/cholmod_l_free_triplet.c",
    "CHOLMOD/Utility/cholmod_transpose.c",
    "CHOLMOD/Utility/cholmod_copy_dense.c",
    "CHOLMOD/Utility/cholmod_l_allocate_dense.c",
    "CHOLMOD/Utility/cholmod_l_allocate_triplet.c",
    "CHOLMOD/Utility/cholmod_l_eye.c",
    "CHOLMOD/Utility/cholmod_copy_dense2.c",
    "CHOLMOD/Utility/cholmod_l_change_factor.c",
    "CHOLMOD/Utility/cholmod_l_pack_factor.c",
    "CHOLMOD/Utility/cholmod_alloc_work.c",
    "CHOLMOD/Utility/cholmod_l_finish.c",
    "CHOLMOD/Utility/cholmod_allocate_sparse.c",
    "CHOLMOD/Utility/cholmod_l_error.c",
    "CHOLMOD/Utility/cholmod_l_realloc_multiple.c",
    "CHOLMOD/Utility/cholmod_l_allocate_work.c",
    "CHOLMOD/Utility/cholmod_l_band_nnz.c",
    "CHOLMOD/Utility/cholmod_speye.c",
    "CHOLMOD/Utility/cholmod_mult_size_t.c",
    "CHOLMOD/Utility/cholmod_l_zeros.c",
    "CHOLMOD/Utility/cholmod_l_free.c",
    "CHOLMOD/Utility/cholmod_cumsum.c",
    "CHOLMOD/Utility/cholmod_l_score_comp.c",
    "CHOLMOD/Utility/cholmod_l_reallocate_sparse.c",
    "CHOLMOD/Utility/cholmod_l_malloc.c",
    "CHOLMOD/Utility/cholmod_l_band.c",
    "CHOLMOD/Utility/cholmod_reallocate_sparse.c",
    "CHOLMOD/Utility/cholmod_sparse_to_dense.c",
    "CHOLMOD/Utility/cholmod_maxrank.c",
    "CHOLMOD/Utility/cholmod_dense_nnz.c",
    "CHOLMOD/Utility/cholmod_transpose_unsym.c",
    "CHOLMOD/Utility/cholmod_realloc.c",
    "CHOLMOD/Utility/cholmod_free_work.c",
    "CHOLMOD/Utility/cholmod_l_add.c",
    "CHOLMOD/Utility/cholmod_copy_factor.c",
    "CHOLMOD/Utility/cholmod_sparse_to_triplet.c",
    "CHOLMOD/Utility/cholmod_l_factor_to_sparse.c",
    "CHOLMOD/Utility/cholmod_band_nnz.c",
    "CHOLMOD/Utility/cholmod_l_divcomplex.c",
    "CHOLMOD/Utility/cholmod_l_sparse_to_triplet.c",
    "CHOLMOD/Utility/cholmod_zeros.c",
    "CHOLMOD/Utility/cholmod_l_triplet_to_sparse.c",
    "CHOLMOD/Utility/cholmod_l_version.c",
    "CHOLMOD/Utility/cholmod_transpose_sym.c",
    "CHOLMOD/Utility/cholmod_l_speye.c",
    "CHOLMOD/Utility/cholmod_ones.c",
    "CHOLMOD/Utility/cholmod_ptranspose.c",
    "CHOLMOD/Utility/cholmod_l_allocate_factor.c",
    "CHOLMOD/Utility/cholmod_free_sparse.c",
    "CHOLMOD/Utility/cholmod_triplet_to_sparse.c",
    "CHOLMOD/Utility/cholmod_l_nnz.c",
    "CHOLMOD/Utility/cholmod_l_copy.c",
    "CHOLMOD/Utility/cholmod_dense_to_sparse.c",
    "CHOLMOD/Utility/cholmod_version.c",
    "CHOLMOD/Utility/cholmod_clear_flag.c",
    "CHOLMOD/Utility/cholmod_l_copy_sparse.c",
    "CHOLMOD/Utility/cholmod_l_set_empty.c",
    "CHOLMOD/Utility/cholmod_l_add_size_t.c",
    "CHOLMOD/Utility/cholmod_error.c",
    "CHOLMOD/Utility/cholmod_l_ensure_dense.c",
    "CHOLMOD/Utility/cholmod_l_sort.c",
    "CHOLMOD/Utility/cholmod_l_realloc.c",
    "CHOLMOD/Utility/cholmod_copy_triplet.c",
    "CHOLMOD/Utility/cholmod_l_maxrank.c",
    "CHOLMOD/Utility/cholmod_aat.c",
    "CHOLMOD/Utility/cholmod_l_dbound.c",
    "CHOLMOD/Utility/cholmod_free_dense.c",
    "CHOLMOD/Utility/cholmod_l_free_factor.c",
    "CHOLMOD/Utility/cholmod_sbound.c",
    "CHOLMOD/Utility/cholmod_memdebug.c",
    "CHOLMOD/Utility/cholmod_l_cumsum.c",
    "CHOLMOD/Utility/cholmod_l_free_dense.c",
    "CHOLMOD/Utility/cholmod_l_spzeros.c",
    "CHOLMOD/Utility/cholmod_free.c",
    "CHOLMOD/Utility/cholmod_copy_sparse.c",
    "CHOLMOD/Utility/cholmod_factor_to_sparse.c",
    "CHOLMOD/Utility/cholmod_add_size_t.c",
    "CHOLMOD/Utility/cholmod_free_factor.c",
    "CHOLMOD/Utility/cholmod_start.c",
    "CHOLMOD/Utility/cholmod_malloc.c",
    "CHOLMOD/Utility/cholmod_l_allocate_sparse.c",
    "CHOLMOD/Utility/cholmod_l_clear_flag.c",
    "CHOLMOD/Utility/cholmod_band.c",
    "CHOLMOD/Utility/cholmod_free_triplet.c",
    "CHOLMOD/Utility/cholmod_l_alloc_factor.c",
    "CHOLMOD/Utility/cholmod_eye.c",
    "CHOLMOD/Utility/cholmod_l_defaults.c",
    "CHOLMOD/Utility/cholmod_calloc.c",
    "CHOLMOD/Utility/cholmod_l_ptranspose.c",
    "CHOLMOD/Utility/cholmod_l_copy_factor.c",
    "CHOLMOD/Utility/cholmod_l_reallocate_triplet.c",
    "CHOLMOD/Utility/cholmod_l_free_sparse.c",
    "CHOLMOD/Utility/cholmod_divcomplex.c",
    "CHOLMOD/Utility/cholmod_spzeros.c",
    "CHOLMOD/Utility/cholmod_l_transpose_sym.c",
    "CHOLMOD/Utility/cholmod_l_transpose.c",
    "CHOLMOD/Utility/cholmod_hypot.c",
    "CHOLMOD/Utility/cholmod_l_sparse_to_dense.c",
    "CHOLMOD/Utility/cholmod_finish.c",
    "CHOLMOD/Utility/cholmod_l_xtype.c",
};

// CHOLMOD flags: build the minimal configuration used by demo_cholmod.c.
const cholmod_flags = [_][]const u8{
    "-DNCHECK",
    "-DNPARTITION",
    "-DNCAMD",
    "-DNMATRIXOPS",
    "-DNMODIFY",
    "-DNSUPERNODAL",
    "-DNPRINT",
};

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const enable_timing = b.option(bool, "timing", "Enable scope timing") orelse false;

    // METIS (+ GKlib), built from vendored source with libc so it compiles on
    // any platform (macOS was previously relying on implicit native libc).
    const gklib_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true });
    gklib_mod.addIncludePath(b.path(metis_dir ++ "/GKlib"));
    gklib_mod.addCSourceFiles(.{
        .root = b.path(metis_dir ++ "/GKlib"),
        .files = &gklib_files,
        // _GNU_SOURCE so glibc declares strptime et al. (macOS declares them
        // unconditionally; glibc gates them behind this feature-test macro).
        .flags = &.{"-D_GNU_SOURCE"},
    });
    const gklib = b.addLibrary(.{ .name = "gk", .linkage = .static, .root_module = gklib_mod });

    const metis_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true });
    metis_mod.addIncludePath(b.path(metis_dir ++ "/include"));
    metis_mod.addIncludePath(b.path(metis_dir ++ "/libmetis"));
    metis_mod.addIncludePath(b.path(metis_dir ++ "/GKlib"));
    metis_mod.addCSourceFiles(.{
        .root = b.path(metis_dir ++ "/libmetis"),
        .files = &metis_files,
        .flags = &.{"-D_GNU_SOURCE"},
    });
    metis_mod.linkLibrary(gklib);
    const libmetis = b.addLibrary(.{ .name = "metis", .linkage = .static, .root_module = metis_mod });

    // Core solver library (C).
    const lib_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true });
    lib_mod.addIncludePath(b.path("src"));
    lib_mod.addIncludePath(b.path(metis_dir ++ "/include"));
    lib_mod.addCSourceFiles(.{
        .files = &.{
            "src/mat.c",
            "src/arena.c",
            "src/linearizer.c",
            "src/solver.c",
            "src/timing.c",
        },
        .flags = &.{},
    });
    const lib = b.addLibrary(.{ .name = "lib", .linkage = .static, .root_module = lib_mod });

    // balTest (C++ reference implementation, uses Eigen).
    const balTest_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libcpp = true });
    balTest_mod.addIncludePath(b.path("src"));
    balTest_mod.addIncludePath(b.path("test/bal"));
    balTest_mod.addIncludePath(b.path(eigen_dir));
    balTest_mod.addCSourceFiles(.{
        .files = &.{
            "test/bal/main.cc",
            "test/bal/sym/rot3.cc",
            "test/bal/sym/ops/rot3/storage_ops.cc",
            "test/bal/sym/ops/rot3/group_ops.cc",
            "test/bal/sym/ops/rot3/lie_group_ops.cc",
            "test/bal/sym/pose3.cc",
            "test/bal/sym/ops/pose3/storage_ops.cc",
            "test/bal/sym/ops/pose3/group_ops.cc",
            "test/bal/sym/ops/pose3/lie_group_ops.cc",
        },
        .flags = &.{},
    });
    balTest_mod.linkLibrary(lib);
    balTest_mod.linkLibrary(libmetis);
    const balTest = b.addExecutable(.{ .name = "balTest", .root_module = balTest_mod });

    // balDemo (the main C demo we optimize / profile).
    const balDemo_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true });
    balDemo_mod.addIncludePath(b.path("src"));
    balDemo_mod.addIncludePath(b.path("test/bal"));
    // -ffast-math doesn't seem to gain much
    // it does turn pow(x, 2) into x * x though
    balDemo_mod.addCSourceFiles(.{
        .files = &.{
            "test/bal/demo.c",
        },
        .flags = &.{},
    });
    balDemo_mod.linkLibrary(lib);
    balDemo_mod.linkLibrary(libmetis);
    if (enable_timing) {
        balDemo_mod.addCMacro("SYM_ENABLE_TIMING", "1");
    }
    const balDemo = b.addExecutable(.{ .name = "balDemo", .root_module = balDemo_mod });

    // CHOLMOD, built from vendored SuiteSparse source.
    const cholmod_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true });
    cholmod_mod.addIncludePath(b.path(suitesparse_dir ++ "/SuiteSparse_config"));
    cholmod_mod.addIncludePath(b.path(suitesparse_dir ++ "/CHOLMOD/Include"));
    cholmod_mod.addCSourceFiles(.{
        .root = b.path(suitesparse_dir),
        .files = &cholmod_files,
        .flags = &cholmod_flags,
    });
    const cholmod = b.addLibrary(.{ .name = "cholmod", .linkage = .static, .root_module = cholmod_mod });

    // balDemoCholmod (alternative solver demo using CHOLMOD).
    const balDemoCholmod_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true });
    balDemoCholmod_mod.addIncludePath(b.path("src"));
    balDemoCholmod_mod.addIncludePath(b.path("test/bal"));
    balDemoCholmod_mod.addIncludePath(b.path(suitesparse_dir ++ "/SuiteSparse_config"));
    balDemoCholmod_mod.addIncludePath(b.path(suitesparse_dir ++ "/CHOLMOD/Include"));
    balDemoCholmod_mod.addCSourceFiles(.{
        .files = &.{
            "test/bal/demo_cholmod.c",
            "test/bal/cholmod_shim.c",
        },
        .flags = &cholmod_flags,
    });
    balDemoCholmod_mod.linkLibrary(lib);
    balDemoCholmod_mod.linkLibrary(libmetis);
    balDemoCholmod_mod.linkLibrary(cholmod);
    if (enable_timing) {
        balDemoCholmod_mod.addCMacro("SYM_ENABLE_TIMING", "1");
    }
    const balDemoCholmod = b.addExecutable(.{ .name = "balDemoCholmod", .root_module = balDemoCholmod_mod });

    // unit tests (C).
    const unit_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true });
    unit_mod.addIncludePath(b.path("src"));
    unit_mod.addCSourceFiles(.{
        .files = &.{
            "test/unit.c",
        },
        .flags = &.{},
    });
    unit_mod.linkLibrary(lib);
    unit_mod.linkLibrary(libmetis);
    const unit = b.addExecutable(.{ .name = "unit", .root_module = unit_mod });

    // Deterministic timing aggregation and scope tests.
    const timingTest_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true });
    timingTest_mod.addIncludePath(b.path("src"));
    timingTest_mod.addCMacro("SYM_ENABLE_TIMING", "1");
    timingTest_mod.addCSourceFiles(.{
        .files = &.{
            "test/timing.c",
        },
        .flags = &.{},
    });
    timingTest_mod.linkLibrary(lib);
    const timingTest = b.addExecutable(.{ .name = "timingTest", .root_module = timingTest_mod });

    b.installArtifact(balTest);
    b.installArtifact(balDemo);
    b.installArtifact(balDemoCholmod);
    b.installArtifact(unit);
    b.installArtifact(timingTest);

    // Named step to build only balDemo (e.g. `zig build balDemo`), avoiding the
    // CHOLMOD/SuiteSparse and Python targets.
    const balDemoStep = b.step("balDemo", "Build only the balDemo executable");
    balDemoStep.dependOn(&b.addInstallArtifact(balDemo, .{}).step);

    // Generate compile_commands.json for clangd with `zig build cdb`.
    const cdb_targets = b.allocator.dupe(*std.Build.Step.Compile, &.{
        gklib,
        libmetis,
        lib,
        balTest,
        balDemo,
        cholmod,
        balDemoCholmod,
        unit,
        timingTest,
    }) catch @panic("OOM");
    _ = zcc.createStep(b, "cdb", cdb_targets);

    // The Python module needs Python + numpy headers at configure time, which
    // aren't available everywhere (e.g. a bare Linux container). Pass -Dnopython
    // to skip it and build only the C/C++ targets.
    const nopython = b.option(bool, "nopython", "Skip the Python module (no Python/numpy needed)") orelse false;
    if (!nopython) {
        const python_exe = b.option([]const u8, "python-exe", "Python executable to use") orelse "python";

        const pythonInc = pythonQuery(b, python_exe, "import sysconfig; print(sysconfig.get_path('include'), end='')");
        const pythonLib = pythonQuery(b, python_exe, "import sysconfig; print(sysconfig.get_config_var('LIBDIR'), end='')");
        const pythonVer = pythonQuery(b, python_exe, "import sysconfig; print(sysconfig.get_config_var('LDVERSION'), end='')");
        const pythonLibName = std.fmt.allocPrint(b.allocator, "python{s}", .{pythonVer}) catch @panic("Missing python");

        const balModule_mod = b.createModule(.{ .target = target, .optimize = optimize, .link_libc = true });
        balModule_mod.addIncludePath(.{ .cwd_relative = pythonInc });
        balModule_mod.addIncludePath(.{
            .cwd_relative = std.fmt.allocPrint(b.allocator, "venv/lib/{s}/site-packages/numpy/core/include/numpy", .{pythonLibName}) catch @panic("Missing python"),
        });
        balModule_mod.addIncludePath(b.path("test/bal"));
        balModule_mod.addIncludePath(b.path("src"));
        balModule_mod.addCSourceFiles(.{
            .files = &.{
                "test/bal/py/balmodule.c",
            },
            .flags = &.{},
        });
        balModule_mod.addLibraryPath(.{ .cwd_relative = pythonLib });
        balModule_mod.linkSystemLibrary(pythonLibName, .{});
        const balModule = b.addLibrary(.{ .name = "balmodule", .linkage = .dynamic, .root_module = balModule_mod });
        // Rename the shared library so Python can find it.
        const balInstallStep = b.addInstallArtifact(balModule, .{ .dest_sub_path = "bal.so" });
        b.getInstallStep().dependOn(&balInstallStep.step);
    }

    const glfw_mod = b.createModule(.{
        .target = target,
        .optimize = optimize,
        .link_libc = true,
    });
    glfw_mod.addCSourceFiles(.{
        .root = b.path("third_party/glfw-3.4/src"),
        .files = &.{
            "cocoa_joystick.m",
            "cocoa_init.m",
            "cocoa_monitor.m",
            "cocoa_time.c",
            "cocoa_window.m",
            "context.c",
            "egl_context.c",
            "glx_context.c",
            "init.c",
            "input.c",
            "monitor.c",
            "nsgl_context.m",
            "null_init.c",
            "null_joystick.c",
            "null_monitor.c",
            "null_window.c",
            "osmesa_context.c",
            "platform.c",
            "posix_module.c",
            "posix_poll.c",
            "posix_thread.c",
            "posix_time.c",
            "vulkan.c",
            "wgl_context.c",
            "window.c",
        },
        .flags = &.{
            "-D_GLFW_COCOA",
        },
    });

    const glfw = b.addLibrary(.{
        .name = "glfw",
        .linkage = .static,
        .root_module = glfw_mod,
    });

    const imgui_mod = b.createModule(.{
        .target = target,
        .optimize = optimize,
        .link_libcpp = true,
    });
    imgui_mod.addIncludePath(b.path("third_party/imgui-1.92.7"));
    imgui_mod.addIncludePath(b.path("third_party/imgui-1.92.7/backends"));
    imgui_mod.addIncludePath(b.path("third_party/glfw-3.4/include"));
    imgui_mod.addCSourceFiles(.{
        .root = b.path("third_party/imgui-1.92.7"),
        .files = &.{
            "imgui.cpp",
            "imgui_demo.cpp",
            "imgui_draw.cpp",
            "imgui_tables.cpp",
            "imgui_widgets.cpp",
            "backends/imgui_impl_glfw.cpp",
            "backends/imgui_impl_metal.mm",
        },
        .flags = &.{},
    });
    imgui_mod.linkLibrary(glfw);

    const imgui = b.addLibrary(.{
        .name = "imgui",
        .linkage = .static,
        .root_module = imgui_mod,
    });

    const main_mod = b.createModule(.{
        .target = target,
        .optimize = optimize,
    });
    main_mod.addIncludePath(b.path("third_party/imgui-1.92.7"));
    main_mod.addIncludePath(b.path("third_party/imgui-1.92.7/backends"));
    main_mod.addIncludePath(b.path("third_party/glfw-3.4/include"));
    main_mod.addCSourceFiles(.{
        .files = &.{
            "src/main.mm",
        },
        .flags = &.{},
    });
    main_mod.linkFramework("Metal", .{});
    main_mod.linkFramework("MetalKit", .{});
    main_mod.linkFramework("Cocoa", .{});
    main_mod.linkFramework("IOKit", .{});
    main_mod.linkFramework("CoreVideo", .{});
    main_mod.linkFramework("QuartzCore", .{});
    main_mod.linkLibrary(imgui);

    const main = b.addExecutable(.{
        .name = "main",
        .root_module = main_mod,
    });

    b.installArtifact(main);
}
