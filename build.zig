const std = @import("std");

pub fn build(b: *std.build.Builder) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // From SymForce:
    // METIS does not have releases recently.  Previously were using nearly the initial commit on
    // github, which is newer than the release on the METIS website. All of the releases on github
    // seem to have memory bugs, which do not appear in this release:
    // https://symforce-org.github.io/downloads/metis-5.1.0.tar.gz
    // SHA256=76faebe03f6c963127dbb73c13eab58c9a3faeae48779f049066a21c087c5db2
    const gklib = b.addStaticLibrary(.{
        .name = "gk",
        .target = target,
        .optimize = optimize,
    });
    gklib.addIncludePath(.{
        .path="metis-5.1.0/GKlib",
    });
    // just `find metis-5.1.0/GKlib -name "*.c"`
    gklib.addCSourceFiles(&.{ 
        "metis-5.1.0/GKlib/htable.c",
        "metis-5.1.0/GKlib/itemsets.c",
        "metis-5.1.0/GKlib/memory.c",
        "metis-5.1.0/GKlib/gkregex.c",
        "metis-5.1.0/GKlib/util.c",
        "metis-5.1.0/GKlib/test/strings.c",
        "metis-5.1.0/GKlib/test/gkgraph.c",
        "metis-5.1.0/GKlib/test/fis.c",
        "metis-5.1.0/GKlib/test/gksort.c",
        "metis-5.1.0/GKlib/test/rw.c",
        "metis-5.1.0/GKlib/omp.c",
        "metis-5.1.0/GKlib/io.c",
        "metis-5.1.0/GKlib/fkvkselect.c",
        "metis-5.1.0/GKlib/sort.c",
        "metis-5.1.0/GKlib/evaluate.c",
        "metis-5.1.0/GKlib/pqueue.c",
        "metis-5.1.0/GKlib/string.c",
        "metis-5.1.0/GKlib/random.c",
        "metis-5.1.0/GKlib/fs.c",
        "metis-5.1.0/GKlib/pdb.c",
        "metis-5.1.0/GKlib/csr.c",
        "metis-5.1.0/GKlib/timers.c",
        "metis-5.1.0/GKlib/error.c",
        "metis-5.1.0/GKlib/seq.c",
        "metis-5.1.0/GKlib/b64.c",
        "metis-5.1.0/GKlib/mcore.c",
        "metis-5.1.0/GKlib/conf/check_thread_storage.c",
        "metis-5.1.0/GKlib/graph.c",
        "metis-5.1.0/GKlib/blas.c",
        "metis-5.1.0/GKlib/getopt.c",
        "metis-5.1.0/GKlib/tokenizer.c",
        "metis-5.1.0/GKlib/rw.c",
    }, &.{});

    const libmetis = b.addStaticLibrary(.{
        .name = "metis",
        .target = target,
        .optimize = optimize,
    });
    libmetis.addIncludePath(.{
        .path="metis-5.1.0/include",
    });
    libmetis.addIncludePath(.{
        .path="metis-5.1.0/libmetis",
    });
    libmetis.addIncludePath(.{
        .path="metis-5.1.0/GKlib",
    });
    libmetis.addCSourceFiles(&.{
        "metis-5.1.0/libmetis/kwayrefine.c",
        "metis-5.1.0/libmetis/mincover.c",
        "metis-5.1.0/libmetis/bucketsort.c",
        "metis-5.1.0/libmetis/parmetis.c",
        "metis-5.1.0/libmetis/util.c",
        "metis-5.1.0/libmetis/kmetis.c",
        "metis-5.1.0/libmetis/meshpart.c",
        "metis-5.1.0/libmetis/compress.c",
        "metis-5.1.0/libmetis/gklib.c",
        "metis-5.1.0/libmetis/auxapi.c",
        "metis-5.1.0/libmetis/separator.c",
        "metis-5.1.0/libmetis/frename.c",
        "metis-5.1.0/libmetis/mcutil.c",
        "metis-5.1.0/libmetis/ometis.c",
        "metis-5.1.0/libmetis/wspace.c",
        "metis-5.1.0/libmetis/sfm.c",
        "metis-5.1.0/libmetis/debug.c",
        "metis-5.1.0/libmetis/balance.c",
        "metis-5.1.0/libmetis/pmetis.c",
        "metis-5.1.0/libmetis/mmd.c",
        "metis-5.1.0/libmetis/refine.c",
        "metis-5.1.0/libmetis/contig.c",
        "metis-5.1.0/libmetis/coarsen.c",
        "metis-5.1.0/libmetis/kwayfm.c",
        "metis-5.1.0/libmetis/stat.c",
        "metis-5.1.0/libmetis/checkgraph.c",
        "metis-5.1.0/libmetis/timing.c",
        "metis-5.1.0/libmetis/fm.c",
        "metis-5.1.0/libmetis/fortran.c",
        "metis-5.1.0/libmetis/initpart.c",
        "metis-5.1.0/libmetis/graph.c",
        "metis-5.1.0/libmetis/srefine.c",
        "metis-5.1.0/libmetis/mesh.c",
        "metis-5.1.0/libmetis/options.c",
        "metis-5.1.0/libmetis/minconn.c",
    }, &.{});
    libmetis.linkLibrary(gklib);

    const lib = b.addStaticLibrary(.{
        .name = "lib",
        .target = target,
        .optimize = optimize,
    });
    lib.addIncludePath(.{
        .path="src",
    });
    lib.addIncludePath(.{
        .path="metis-5.1.0/include",
    });
    lib.addCSourceFiles(&.{
        "src/mat.c",
        "src/arena.c",
        "src/linearizer.c",
        "src/solver.c",
    }, &.{});

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
    balTest.addIncludePath(.{
        .path="eigen-3.4.0/include/eigen3",
    });
    balTest.addCSourceFiles(&.{
        "test/bal/main.cc",
        "test/bal/sym/rot3.cc",
        "test/bal/sym/ops/rot3/storage_ops.cc",
        "test/bal/sym/ops/rot3/group_ops.cc",
        "test/bal/sym/ops/rot3/lie_group_ops.cc",
        "test/bal/sym/pose3.cc",
        "test/bal/sym/ops/pose3/storage_ops.cc",
        "test/bal/sym/ops/pose3/group_ops.cc",
        "test/bal/sym/ops/pose3/lie_group_ops.cc",
    }, &.{
    });
    balTest.linkLibrary(lib);
    balTest.linkLibrary(libmetis);
    balTest.linkLibCpp();

    const unit = b.addExecutable(.{
        .name = "unit",
        .target = target,
        .optimize = optimize,
    });
    unit.addIncludePath(.{
        .path="src",
    });
    unit.addCSourceFiles(&.{
        "test/unit.c",
    }, &.{});
    unit.linkLibrary(lib);
    unit.linkLibrary(libmetis);

    b.installArtifact(balTest);
    b.installArtifact(unit);
}
