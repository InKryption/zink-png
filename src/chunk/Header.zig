const std = @import("std");

const Header = @This();
length: u31 align(@alignOf(u32)),
type: Type align(@alignOf(u32)),

pub const Type = enum(u32) {
    pub const Tag = @typeInfo(Type).Enum.tag_type;
    // Standard Critical Chunk Types
    IHDR = std.mem.readIntBig(u32, "IHDR"),
    PLTE = std.mem.readIntBig(u32, "PLTE"),
    IDAT = std.mem.readIntBig(u32, "IDAT"),
    IEND = std.mem.readIntBig(u32, "IEND"),

    // Standard Ancillary Chunk Types
    bKGD = std.mem.readIntBig(u32, "bKGD"),
    cHRM = std.mem.readIntBig(u32, "cHRM"),
    gAMA = std.mem.readIntBig(u32, "gAMA"),
    hIST = std.mem.readIntBig(u32, "hIST"),
    pHYs = std.mem.readIntBig(u32, "pHYs"),
    sBIT = std.mem.readIntBig(u32, "sBIT"),
    tEXt = std.mem.readIntBig(u32, "tEXt"),
    tIME = std.mem.readIntBig(u32, "tIME"),
    tRNS = std.mem.readIntBig(u32, "tRNS"),
    zTXt = std.mem.readIntBig(u32, "zTXt"),

    // Non-Standard Chunk Types
    _,

    pub fn from(bytes: [4]u8) Type {
        return @intToEnum(Type, (std.mem.readIntBig(u32, bytes[0..])));
    }

    pub fn string(self: Type) [4]u8 {
        return std.mem.toBytes(self.intBig());
    }

    pub fn isValidAscii(self: Type) bool {
        for (self.string()) |byte| {
            if (!std.ascii.isAlpha(byte)) return false;
        }
        return true;
    }

    pub fn int(self: Type) Tag {
        return @enumToInt(self);
    }
    pub fn intBig(self: Type) Tag {
        return std.mem.nativeToBig(Tag, self.int());
    }
    pub fn intLittle(self: Type) Tag {
        return std.mem.nativeToLittle(Tag, self.int());
    }

    pub fn property(self: Type, byte_index: u2) bool {
        return (self.string()[byte_index] & 32) != 0;
    }

    pub fn format(
        ch_type: Type,
        comptime fmt_str: []const u8,
        fmt_options: std.fmt.FormatOptions,
        writer: anytype,
    ) @TypeOf(writer).Error!void {
        _ = fmt_options;
        const lazy = struct {
            const unknown_format_str = @compileError("Unknown format string '" ++ fmt_str ++ "'.\n");
        };

        switch (fmt_str.len) {
            0 => {},
            else => lazy.unknown_format_str,
        }
        return writer.print("{s}({s})", .{ @typeName(@This()), ch_type.string() });
    }
};

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
    pub const Type = union(enum) { ok: Header.Type, invalid: Header.Type };

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
            const ch_type = Type.from(fbs.reader().readBytesNoEof(4) catch unreachable);
            break :ch_type if (ch_type.isValidAscii())
                FromBytesResult.Type{ .ok = ch_type }
            else
                FromBytesResult.Type{ .invalid = ch_type };
        },
    };
}

test "Type & Header" {
    try std.testing.expectEqualStrings("WOAH", &Type.from("WOAH".*).string());
    for (@as([4096]void, undefined)) |_, sample_valid_ch_len| {
        inline for (comptime std.enums.values(Type)) |sample_valid_ch_ty| {
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
                .type = .{ .invalid = Type.from("0000".*) },
            },
            Header.fromBytes(blk: {
                var bytes_with_invalid_type: [8]u8 = undefined;
                bytes_with_invalid_type[0..4].* = std.mem.toBytes(std.mem.nativeToBig(u32, @intCast(u31, sample_valid_ch_len)));
                bytes_with_invalid_type[4..8].* = Type.from("0000".*).string();
                break :blk bytes_with_invalid_type;
            }),
        );
    }
}
