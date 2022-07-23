const std = @import("std");
const png = @import("../main.zig");

const Header = @This();
length: u31 align(@alignOf(u32)),
type: png.chunk.Type align(@alignOf(u32)),

pub fn toBytes(header: Header) [8]u8 {
    var result: [8]u8 = undefined;
    var fbs = std.io.fixedBufferStream(result[0..]);
    fbs.writer().writeIntBig(u32, header.length) catch unreachable;
    fbs.writer().writeAll(header.type.string()[0..]) catch unreachable;
    return result;
}

pub const FromBytesResult = struct {
    length: Length,
    type: FromBytesResult.Type,

    pub const Length = union(enum) { ok: u31, invalid: u32 };
    /// `ok.isValidAscii()` => `true`
    /// `invalid.isValidAscii()` => `false`
    pub const Type = union(enum) { ok: png.chunk.Type, invalid: png.chunk.Type };

    pub const UnwrapError = error{
        InvalidLength,
        InvalidType,
    };
    pub fn unwrap(result: FromBytesResult) UnwrapError!Header {
        return Header{
            .length = switch (result.length) {
                .ok => |length| length,
                .invalid => return error.InvalidLength,
            },
            .type = switch (result.type) {
                .ok => |ch_type| ch_type,
                .invalid => return error.InvalidType,
            },
        };
    }
};

pub fn fromBytes(bytes: [8]u8) FromBytesResult {
    var fbs = std.io.fixedBufferStream(bytes[0..]);
    return FromBytesResult{
        .length = length: {
            const raw_len = fbs.reader().readIntBig(u32) catch unreachable;
            break :length if (std.math.cast(u31, raw_len)) |length|
                FromBytesResult.Length{ .ok = length }
            else
                FromBytesResult.Length{ .invalid = raw_len };
        },
        .type = ch_type: {
            const ch_type = png.chunk.Type.from(fbs.reader().readBytesNoEof(4) catch unreachable);
            break :ch_type if (ch_type.isValidAscii())
                FromBytesResult.Type{ .ok = ch_type }
            else
                FromBytesResult.Type{ .invalid = ch_type };
        },
    };
}

test "Type & Header" {
    try std.testing.expectEqualStrings("WOAH", &png.chunk.Type.from("WOAH".*).string());
    for (@as([4096]void, undefined)) |_, sample_valid_ch_len| {
        inline for (comptime std.enums.values(png.chunk.Type)) |sample_valid_ch_ty| {
            const valid_header = Header{
                .length = @intCast(u31, sample_valid_ch_len),
                .type = sample_valid_ch_ty,
            };
            try std.testing.expectEqual(valid_header, Header.fromBytes(valid_header.toBytes()).unwrap() catch @panic("That shouldn't have happened."));

            const invalid_len: u32 = std.math.maxInt(u31) + 1 + @intCast(u32, sample_valid_ch_len);
            try std.testing.expectEqual(Header.FromBytesResult{
                .length = .{ .invalid = invalid_len },
                .type = .{ .ok = sample_valid_ch_ty },
            }, Header.fromBytes(blk: {
                var bytes_with_invalid_len: [8]u8 = undefined;
                bytes_with_invalid_len[0..4].* = std.mem.toBytes(std.mem.nativeToBig(u32, invalid_len));
                bytes_with_invalid_len[4..8].* = sample_valid_ch_ty.string();
                break :blk bytes_with_invalid_len;
            }));
        }
        try std.testing.expectEqual(
            Header.FromBytesResult{
                .length = .{ .ok = @intCast(u31, sample_valid_ch_len) },
                .type = .{ .invalid = png.chunk.Type.from("0000".*) },
            },
            Header.fromBytes(blk: {
                var bytes_with_invalid_type: [8]u8 = undefined;
                bytes_with_invalid_type[0..4].* = std.mem.toBytes(std.mem.nativeToBig(u32, @intCast(u31, sample_valid_ch_len)));
                bytes_with_invalid_type[4..8].* = png.chunk.Type.from("0000".*).string();
                break :blk bytes_with_invalid_type;
            }),
        );
    }
}
