const std = @import("std");

pub fn getPackage(b: *std.build.Builder, name: []const u8) std.build.Pkg {
    return std.build.Pkg{
        .name = name,
        .source = .{ .path = b.pathJoin(std.fs.path.dirname(@src().file).?, "src/main.zig") },
        .dependencies = null,
    };
}

pub fn build(b: *std.build.Builder) void {
    const mode = b.standardReleaseOptions();

    const main_tests = b.addTest("src/main.zig");
    main_tests.setBuildMode(mode);

    const test_step = b.step("test", "Run library tests");
    test_step.dependOn(&main_tests.step);
}
