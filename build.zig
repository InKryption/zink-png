const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const lib_module = b.createModule(.{ .source_file = .{ .path = "src/main.zig" } });
    b.addModule(.{
        .name = "zink-png",
        .source_file = lib_module.source_file,
    });

    const unit_tests = b.addTest(.{
        .root_source_file = lib_module.source_file,
        .target = target,
        .optimize = optimize,
    });
    const unit_tests_step = b.step("unit-tests", "Run unit tests");
    unit_tests_step.dependOn(&unit_tests.step);

    const test_cases = b.addTest(.{
        .root_source_file = .{ .path = "test/cases.zig" },
        .target = target,
        .optimize = optimize,
    });
    test_cases.addModule("zink-png", lib_module);

    const test_cases_step = b.step("test-cases", "Run test cases");
    test_cases_step.dependOn(&test_cases.step);
}
