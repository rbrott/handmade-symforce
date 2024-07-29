const builtin = @import("builtin");
const std = @import("std");

// cribbed from https://github.com/spiraldb/ziggy-pydust/blob/16da6a0cc4ec10295526b7abef0bcfb9dabb65f0/build.zig
const runProcess = if (builtin.zig_version.minor >= 12) std.process.Child.run else std.process.Child.exec;

fn getPythonIncludePath(
    python_exe: []const u8,
    allocator: std.mem.Allocator,
) ![]const u8 {
    const includeResult = try runProcess(.{
        .allocator = allocator,
        .argv = &.{ python_exe, "-c", "import sysconfig; print(sysconfig.get_path('include'), end='')" },
    });
    defer allocator.free(includeResult.stderr);
    return includeResult.stdout;
}

fn getPythonLibraryPath(python_exe: []const u8, allocator: std.mem.Allocator) ![]const u8 {
    const includeResult = try runProcess(.{
        .allocator = allocator,
        .argv = &.{ python_exe, "-c", "import sysconfig; print(sysconfig.get_config_var('LIBDIR'), end='')" },
    });
    defer allocator.free(includeResult.stderr);
    return includeResult.stdout;
}

fn getPythonLDVersion(python_exe: []const u8, allocator: std.mem.Allocator) ![]const u8 {
    const includeResult = try runProcess(.{
        .allocator = allocator,
        .argv = &.{ python_exe, "-c", "import sysconfig; print(sysconfig.get_config_var('LDVERSION'), end='')" },
    });
    defer allocator.free(includeResult.stderr);
    return includeResult.stdout;
}

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const metis = b.dependency("metis", .{
        .target = target,
        .optimize = optimize,
    });
    const libmetis = metis.artifact("metis");

    const lib = b.addStaticLibrary(.{
        .name = "lib",
        .target = target,
        .optimize = optimize,
    });
    lib.addIncludePath(b.path("src"));
    lib.addIncludePath(metis.path("include"));
    lib.addCSourceFiles(.{
        .files = &.{
            "src/mat.c",
            "src/arena.c",
            "src/linearizer.c",
            "src/solver.c",
        },
        .flags = &.{}
    });

    const eigen = b.dependency("eigen", .{});

    const balTest = b.addExecutable(.{
        .name = "balTest",
        .target = target,
        .optimize = optimize,
    });
    balTest.addIncludePath(b.path("src"));
    balTest.addIncludePath(b.path("test/bal"));
    balTest.addIncludePath(eigen.path(""));
    balTest.addCSourceFiles(.{
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
        .flags = &.{}
    });
    balTest.linkLibrary(lib);
    balTest.linkLibrary(libmetis);
    balTest.linkLibCpp();

    const balDemo = b.addExecutable(.{
        .name = "balDemo",
        .target = target,
        .optimize = optimize,
    });
    balDemo.addIncludePath(b.path("src"));
    balDemo.addIncludePath(b.path("test/bal"));
    // -ffast-math doesn't seem to gain much
    // it does turn pow(x, 2) into x * x though
    balDemo.addCSourceFiles(.{
        .files = &.{
            "test/bal/demo.c",
        },
        .flags = &.{}
    });
    balDemo.linkLibrary(lib);
    balDemo.linkLibrary(libmetis);

    const suitesparse = b.dependency("suitesparse", .{
        .target = target,
        .optimize = optimize,
    });
    const cholmod = suitesparse.artifact("cholmod");

    const balDemoCholmod = b.addExecutable(.{
        .name = "balDemoCholmod",
        .target = target,
        .optimize = optimize,
    });
    balDemoCholmod.addIncludePath(b.path("src"));
    balDemoCholmod.addIncludePath(b.path("test/bal"));
    balDemoCholmod.addIncludePath(suitesparse.path("SuiteSparse_config"));
    balDemoCholmod.addIncludePath(suitesparse.path("CHOLMOD/Include"));
    balDemoCholmod.addCSourceFiles(.{
        .files = &.{
            "test/bal/demo_cholmod.c",
            "test/bal/cholmod_shim.c",
        },
        .flags = &.{
            "-DNCHECK",
            "-DNPARTITION",
            "-DNCAMD",
            "-DNMATRIXOPS",
            "-DNMODIFY",
            "-DNSUPERNODAL",
            "-DNPRINT",
        },
    });
    balDemoCholmod.linkLibrary(lib);
    balDemoCholmod.linkLibrary(libmetis);
    balDemoCholmod.linkLibrary(cholmod);

    const unit = b.addExecutable(.{
        .name = "unit",
        .target = target,
        .optimize = optimize,
    });
    unit.addIncludePath(b.path("src"));
    unit.addCSourceFiles(.{
        .files = &.{
            "test/unit.c",
        },
        .flags = &.{}
    });
    unit.linkLibrary(lib);
    unit.linkLibrary(libmetis);

    b.installArtifact(balTest);
    b.installArtifact(balDemo);
    b.installArtifact(balDemoCholmod);
    b.installArtifact(unit);

    const python_exe = b.option([]const u8, "python-exe", "Python executable to use") orelse "python";

    const pythonInc = getPythonIncludePath(python_exe, b.allocator) catch @panic("Missing python");
    const pythonLib = getPythonLibraryPath(python_exe, b.allocator) catch @panic("Missing python");
    const pythonVer = getPythonLDVersion(python_exe, b.allocator) catch @panic("Missing python");
    const pythonLibName = std.fmt.allocPrint(b.allocator, "python{s}", .{pythonVer}) catch @panic("Missing python");

    const balModule = b.addSharedLibrary(.{
        .name = "balmodule",
        .target = target,
        .optimize = optimize,
    });
    balModule.addIncludePath(.{ .cwd_relative = pythonInc });
    balModule.addIncludePath(.{
        .cwd_relative = std.fmt.allocPrint(b.allocator, "venv/lib/{s}/site-packages/numpy/core/include/numpy", .{pythonLibName}) catch @panic("Missing python"),
    });
    balModule.addIncludePath(b.path("test/bal"));
    balModule.addIncludePath(b.path("src"));
    balModule.addCSourceFiles(.{
        .files = &.{
            "test/bal/py/balmodule.c",
        },
        .flags = &.{}
    });
    balModule.addLibraryPath(.{ .cwd_relative = pythonLib });
    balModule.linkSystemLibrary(pythonLibName);
    // Rename the shared library so Python can find it.
    const balInstallStep = b.addInstallArtifact(balModule, .{ .dest_sub_path = "bal.so" });
    b.getInstallStep().dependOn(&balInstallStep.step);
}
