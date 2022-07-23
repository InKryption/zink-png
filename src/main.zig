const std = @import("std");
const builtin = @import("builtin");

const util = @import("util.zig");

pub const signature: [8]u8 = [8]u8{ 137, 80, 78, 71, 13, 10, 26, 10 };
pub const chunk = @import("chunk/chunk.zig");

comptime {
    _ = chunk;
}
