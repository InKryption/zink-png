const std = @import("std");

const PLTE = @This();
entries_buf: [256]Entry,
count_minus_one: u8,

pub const Entry = extern struct { r: u8, g: u8, b: u8 };
pub fn count(plte: PLTE) std.math.IntFittingRange(1, 256) {
    return @as(std.math.IntFittingRange(1, 256), plte.count_minus_one) + 1;
}
pub fn entries(plte: *const PLTE) []const Entry {
    return plte.entries_buf[0..plte.count()];
}

pub const RawBytes = std.BoundedArray(u8, 256 * 3);

pub fn toBytes(plte: PLTE) RawBytes {
    var result = RawBytes{};
    for (plte.entries()) |entry| {
        result.appendSlice(&[_]u8{ entry.r, entry.g, entry.b }) catch unreachable;
    }
    return result;
}

pub const FromBytesResult = union(enum) {
    ok: PLTE,
    missing_one_byte: MissingOneByte,
    missing_two_bytes: MissingTwoBytes,

    pub const MissingOneByte = struct {
        full: std.BoundedArray(Entry, 255),
        trailing: Trailing,
        pub const Trailing = struct { r: u8, g: u8 };
    };
    pub const MissingTwoBytes = struct {
        full: std.BoundedArray(Entry, 255),
        trailing: Trailing,
        pub const Trailing = struct { r: u8 };
    };

    pub const UnwrapError = error{
        MissingOneByte,
        MissingTwoBytes,
    };
    pub fn unwrap(result: FromBytesResult) UnwrapError!PLTE {
        return switch (result) {
            .ok => |value| value,
            .missing_one_byte => error.MissingOneByte,
            .missing_two_bytes => error.MissingTwoBytes,
        };
    }
};
/// Asserts `bytes.len > 0`.
pub fn fromBytes(bytes: RawBytes) FromBytesResult {
    std.debug.assert(bytes.len > 0);

    var result_entries = std.BoundedArray(Entry, 256){};
    const used_bytes_len = bytes.len - (bytes.len % 3);

    comptime std.debug.assert(@typeInfo(Entry).Struct.layout == .Extern);
    comptime std.debug.assert(@typeInfo(Entry).Struct.fields.len == 3);
    comptime std.debug.assert(std.mem.eql(u8, @typeInfo(Entry).Struct.fields[0].name, "r"));
    comptime std.debug.assert(std.mem.eql(u8, @typeInfo(Entry).Struct.fields[1].name, "g"));
    comptime std.debug.assert(std.mem.eql(u8, @typeInfo(Entry).Struct.fields[2].name, "b"));
    result_entries.appendSliceAssumeCapacity(std.mem.bytesAsSlice(Entry, bytes.constSlice()[0..used_bytes_len]));

    return switch (bytes.len - used_bytes_len) {
        0 => FromBytesResult{ .ok = PLTE{
            .entries_buf = result_entries.buffer,
            .count_minus_one = @intCast(u8, result_entries.len - 1),
        } },
        1 => FromBytesResult{ .missing_one_byte = FromBytesResult.MissingOneByte{
            .full = std.BoundedArray(Entry, 255).fromSlice(result_entries.constSlice()) catch unreachable,
            .trailing = .{
                .r = bytes.get(bytes.len - 2),
                .g = bytes.get(bytes.len - 1),
            },
        } },
        2 => FromBytesResult{ .missing_two_bytes = FromBytesResult.MissingTwoBytes{
            .full = std.BoundedArray(Entry, 255).fromSlice(result_entries.constSlice()) catch unreachable,
            .trailing = .{
                .r = bytes.get(bytes.len - 1),
            },
        } },
        else => unreachable,
    };
}

test "PLTE" {
    var plte = comptime PLTE.fromBytes(PLTE.RawBytes.fromSlice(std.mem.asBytes(&[_][3]u8{
        //  r    g    b
        .{ 000, 000, 000 },
        .{ 255, 000, 000 },
        .{ 000, 255, 000 },
        .{ 000, 000, 255 },
        .{ 000, 255, 255 },
        .{ 255, 000, 255 },
        .{ 255, 255, 000 },
        .{ 255, 255, 255 },
    })) catch @panic("Oof, didn't meant to do that")).ok;
    try std.testing.expectEqual(plte, PLTE.fromBytes(plte.toBytes()).unwrap() catch @panic("Something went wrong there."));
}
