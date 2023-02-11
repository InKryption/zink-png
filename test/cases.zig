const std = @import("std");
const builtin = @import("builtin");

pub const png = @import("zink-png");

comptime {
    std.debug.assert(builtin.is_test);
    _ = @import("cases/basic.zig");
}
