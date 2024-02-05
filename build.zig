const std = @import("std");

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
    lib.addIncludePath(.{
        .path="src",
    });
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
    balTest.addIncludePath(.{
        .path="src",
    });
    balTest.addIncludePath(.{
        .path="test/bal",
    });
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
    balDemo.addIncludePath(.{
        .path="src",
    });
    balDemo.addIncludePath(.{
        .path="test/bal",
    });
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
    balDemoCholmod.addIncludePath(.{
        .path="src",
    });
    balDemoCholmod.addIncludePath(.{
        .path="test/bal",
    });
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
    unit.addIncludePath(.{
        .path="src",
    });
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
}
