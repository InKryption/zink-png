const builtin = @import("builtin");

pub const signature = [_]u8{ 137, 80, 78, 71, 13, 10, 26, 10 };
pub const chunk = @import("chunk.zig");

comptime {
    if (builtin.is_test) {
        _ = chunk;
    }
}
