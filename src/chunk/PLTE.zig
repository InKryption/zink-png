const std = @import("std");
const chunk = @import("chunk.zig");

const PLTE = @This();
entries_buf: [256]Entry,
count_minus_one: u8,

/// r: entry[0]
/// g: entry[1]
/// b: entry[2]
pub const Entry = [3]u8;
pub fn count(plte: PLTE) std.math.IntFittingRange(1, 256) {
    return @as(std.math.IntFittingRange(1, 256), plte.count_minus_one) + 1;
}
pub fn entries(plte: *const PLTE) []const Entry {
    return plte.entries_buf[0..plte.count()];
}

pub const RawBytes = std.BoundedArray(u8, 256 * 3);

pub fn toBytes(plte: PLTE) RawBytes {
    var result = RawBytes{};
    result.appendSliceAssumeCapacity(std.mem.sliceAsBytes(plte.entries()));
    return result;
}

pub const FromBytesResult = union(enum) {
    ok: PLTE,
    missing_one_byte: MissingOneByte,
    missing_two_bytes: MissingTwoBytes,

    pub const MissingOneByte = struct {
        full: Full,
        trailing: Trailing,
        pub const Full = std.BoundedArray(Entry, 255);
        pub const Trailing = struct { r: u8, g: u8 };
    };
    pub const MissingTwoBytes = struct {
        full: Full,
        trailing: Trailing,
        pub const Full = std.BoundedArray(Entry, 255);
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

    result_entries.appendSliceAssumeCapacity(std.mem.bytesAsSlice(Entry, bytes.constSlice()[0..used_bytes_len]));

    return switch (bytes.len - used_bytes_len) {
        0 => FromBytesResult{ .ok = PLTE{
            .entries_buf = result_entries.buffer,
            .count_minus_one = @intCast(u8, result_entries.len - 1),
        } },
        1 => FromBytesResult{ .missing_two_bytes = FromBytesResult.MissingTwoBytes{
            .full = std.BoundedArray(Entry, 255).fromSlice(result_entries.constSlice()) catch unreachable,
            .trailing = .{
                .r = bytes.get(bytes.len - 1),
            },
        } },
        2 => FromBytesResult{ .missing_one_byte = FromBytesResult.MissingOneByte{
            .full = std.BoundedArray(Entry, 255).fromSlice(result_entries.constSlice()) catch unreachable,
            .trailing = .{
                .r = bytes.get(bytes.len - 2),
                .g = bytes.get(bytes.len - 1),
            },
        } },
        else => unreachable,
    };
}

pub const chunk_type: chunk.Header.Type = .PLTE;
pub fn calculateCrc(ihdr: PLTE) u32 {
    var hasher = std.hash.Crc32.init();
    hasher.hash(&chunk.Header.Type.string(chunk_type));
    hasher.hash(&ihdr.toBytes());
    return hasher.final();
}

test "PLTE" {
    const correct_bytes = @bitCast([8 * 3]u8, [_][3]u8{
        //  r    g    b
        .{ 000, 000, 000 },
        .{ 255, 000, 000 },

        .{ 000, 255, 000 },
        .{ 000, 000, 255 },

        .{ 000, 255, 255 },
        .{ 255, 000, 255 },

        .{ 255, 255, 000 },
        .{ 255, 255, 255 },
    });

    var plte = PLTE.fromBytes(PLTE.RawBytes.fromSlice(correct_bytes[0..]) catch @panic("Oof, didn't meant to do that")).ok;
    try std.testing.expectEqual(plte, PLTE.fromBytes(plte.toBytes()).unwrap() catch @panic("Something went wrong there."));

    try std.testing.expectEqual(PLTE.FromBytesResult{
        .missing_one_byte = PLTE.FromBytesResult.MissingOneByte{
            .full = PLTE.FromBytesResult.MissingOneByte.Full.fromSlice(plte.entries()) catch unreachable,
            .trailing = .{ .r = 42, .g = 7 },
        },
    }, PLTE.fromBytes(PLTE.RawBytes.fromSlice(&correct_bytes ++ [_]u8{ 42, 7 }) catch unreachable));

    try std.testing.expectEqual(PLTE.FromBytesResult{
        .missing_two_bytes = PLTE.FromBytesResult.MissingTwoBytes{
            .full = PLTE.FromBytesResult.MissingTwoBytes.Full.fromSlice(plte.entries()) catch unreachable,
            .trailing = .{ .r = 42 },
        },
    }, PLTE.fromBytes(PLTE.RawBytes.fromSlice(&correct_bytes ++ [_]u8{42}) catch unreachable));
}
