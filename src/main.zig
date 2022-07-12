const std = @import("std");
const builtin = @import("builtin");

const util = @import("util.zig");

pub const signature: [8]u8 = .{ 137, 80, 78, 71, 13, 10, 26, 10 };

pub const ChunkType = enum(u32) {
    pub const Tag = @typeInfo(ChunkType).Enum.tag_type;
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

    pub fn from(bytes: *const [4]u8) ChunkType {
        return @intToEnum(ChunkType, (std.mem.readIntBig(u32, bytes)));
    }

    pub fn string(self: ChunkType) [4]u8 {
        return std.mem.toBytes(self.intBig());
    }

    pub fn isValidAscii(self: ChunkType) bool {
        for (self.string()) |byte| {
            if (!std.ascii.isAlpha(byte)) return false;
        }
        return true;
    }

    pub fn int(self: ChunkType) Tag {
        return @enumToInt(self);
    }
    pub fn intBig(self: ChunkType) Tag {
        return std.mem.nativeToBig(Tag, self.int());
    }
    pub fn intLittle(self: ChunkType) Tag {
        return std.mem.nativeToLittle(Tag, self.int());
    }

    pub fn property(self: ChunkType, byte_index: u2) bool {
        return (self.string()[byte_index] & 32) != 0;
    }

    pub fn format(
        ch_type: ChunkType,
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

pub const ChunkHeader = struct {
    length: u31 align(@alignOf(u32)),
    type: ChunkType align(@alignOf(u32)),

    pub fn toBytes(header: ChunkHeader) [8]u8 {
        var result: [8]u8 = undefined;
        result[0..4].* = std.mem.toBytes(std.mem.nativeToBig(u32, header.length));
        result[4..8].* = header.type.string();
        return result;
    }

    pub const FromBytesResult = struct {
        length: Length,
        type: ChunkType,

        pub const Length = union(enum) {
            ok: u31,
            invalid: u32,
        };

        pub const UnwrapError = error{InvalidLength};
        pub fn unwrap(result: FromBytesResult) UnwrapError!ChunkHeader {
            return ChunkHeader{
                .length = switch (result.length) {
                    .ok => |length| length,
                    .invalid => return error.InvalidLength,
                },
                .type = result.type,
            };
        }
    };

    pub fn fromBytes(bytes: [8]u8) FromBytesResult {
        var fbs = std.io.fixedBufferStream(bytes[0..]);
        var result: FromBytesResult = undefined;

        result.length = length: {
            const raw_len = fbs.reader().readIntBig(u32) catch unreachable;
            break :length if (std.math.cast(u31, raw_len)) |length|
                FromBytesResult.Length{ .ok = length }
            else
                FromBytesResult.Length{ .invalid = raw_len };
        };

        result.type = ChunkType.from(&(fbs.reader().readBytesNoEof(4) catch unreachable));

        std.debug.assert(fbs.pos == bytes.len);
        return result;
    }
};

test "ChunkType & ChunkHeader" {
    try std.testing.expectEqualStrings("WOAH", &ChunkType.from("WOAH").string());
    for (@as([4096]void, undefined)) |_, i| {
        inline for (comptime std.enums.values(ChunkType)) |ch_type| {
            const valid_header = ChunkHeader{
                .length = @intCast(u31, i),
                .type = ch_type,
            };
            try std.testing.expectEqual(valid_header, ChunkHeader.fromBytes(valid_header.toBytes()).unwrap() catch @panic("That shouldn't have happened."));

            const invalid_len: u32 = comptime std.math.maxInt(u31) + 1 + @intCast(u32, i);
            const bytes_with_invalid_len: [8]u8 = blk: {
                var bytes_with_invalid_len: [8]u8 = undefined;
                bytes_with_invalid_len[0..4].* = std.mem.toBytes(std.mem.nativeToBig(u32, invalid_len));
                bytes_with_invalid_len[4..8].* = ch_type.string();
                break :blk bytes_with_invalid_len;
            };
            try std.testing.expectEqual(ChunkHeader.FromBytesResult{
                .length = .{ .invalid = invalid_len },
                .type = ch_type,
            }, ChunkHeader.fromBytes(bytes_with_invalid_len));
        }
    }
}

pub const ChunkDataIHDR = struct {
    width: u31 align(@alignOf(u32)),
    height: u31 align(@alignOf(u32)),

    bit_depth: BitDepth align(@alignOf(u8)),
    color_type: ColorType align(@alignOf(u8)),

    compression_method: CompressionMethod align(@alignOf(u8)),
    filter_method: FilterMethod align(@alignOf(u8)),
    interlace_method: InterlaceMethod align(@alignOf(u8)),

    pub const BitDepth = enum(u4) {
        @"1" = 1,
        @"2" = 2,
        @"4" = 4,
        @"8" = 8,
        @"16" = 16,
    };
    pub const ColorType = enum(u3) {
        // zig fmt: off
        grayscale          = 0 + 0 + 0, // 0
        rgb         = 0 + 2 + 0, // 2
        palette = 1 + 2 + 0, // 3
        grayscale_alpha         = 0 + 0 + 4, // 4
        rgb_alpha   = 0 + 2 + 4, // 6
        // zig fmt: on

        pub fn colorEnabled(color_type: ColorType) bool {
            return @enumToInt(color_type) & 2 != 0;
        }
        pub fn alphaEnabled(color_type: ColorType) bool {
            return @enumToInt(color_type) & 4 != 0;
        }
        pub fn paletteEnabled(color_type: ColorType) bool {
            return @enumToInt(color_type) & 1 != 0;
        }

        pub fn bitDepthEnabled(color_type: ColorType, bit_depth: BitDepth) bool {
            return switch (color_type) {
                .grayscale => switch (bit_depth) {
                    .@"1",
                    .@"2",
                    .@"4",
                    .@"8",
                    .@"16",
                    => true,
                },
                .rgb => switch (bit_depth) {
                    .@"1",
                    .@"2",
                    .@"4",
                    => false,
                    .@"8",
                    .@"16",
                    => true,
                },
                .palette => switch (bit_depth) {
                    .@"1",
                    .@"2",
                    .@"4",
                    .@"8",
                    => true,
                    .@"16" => false,
                },
                .grayscale_alpha => switch (bit_depth) {
                    .@"1",
                    .@"2",
                    .@"4",
                    => false,
                    .@"8",
                    .@"16",
                    => true,
                },
                .rgb_alpha => switch (bit_depth) {
                    .@"1",
                    .@"2",
                    .@"4",
                    => false,
                    .@"8",
                    .@"16",
                    => true,
                },
            };
        }
    };

    pub const CompressionMethod = enum(u8) {
        @"0" = 0,
        _,
    };
    pub const FilterMethod = enum(u8) {
        @"0" = 0,
        _,
    };
    pub const InterlaceMethod = enum(u8) {
        none = 0,
        adam7 = 1,
        _,
    };

    pub fn toBytes(ihdr: ChunkDataIHDR) [13]u8 {
        return [0]u8{} ++
            std.mem.toBytes(std.mem.nativeToBig(u32, ihdr.width)) ++
            std.mem.toBytes(std.mem.nativeToBig(u32, ihdr.height)) ++
            std.mem.toBytes(std.mem.nativeToBig(u8, @enumToInt(ihdr.bit_depth))) ++
            std.mem.toBytes(std.mem.nativeToBig(u8, @enumToInt(ihdr.color_type))) ++
            std.mem.toBytes(std.mem.nativeToBig(u8, @enumToInt(ihdr.compression_method))) ++
            std.mem.toBytes(std.mem.nativeToBig(u8, @enumToInt(ihdr.filter_method))) ++
            std.mem.toBytes(std.mem.nativeToBig(u8, @enumToInt(ihdr.interlace_method)));
    }

    pub const FromBytesResult = struct {
        width: WidthResult,
        height: HeightResult,

        bit_depth: BitDepthResult,
        color_type: ColorTypeResult,

        compression_method: CompressionMethod,
        filter_method: FilterMethod,
        interlace_method: InterlaceMethod,

        pub const WidthResult = union(enum) { ok: u31, invalid: u32 };
        pub const HeightResult = union(enum) { ok: u31, invalid: u32 };

        pub const BitDepthResult = union(enum) { ok: BitDepth, invalid: u8, invalid_for_color: BitDepth };
        pub const ColorTypeResult = union(enum) { ok: ColorType, invalid: u8 };

        pub const UnwrapError = error{
            InvalidWidth,
            InvalidHeight,
            InvalidBitDepth,
            InvalidColorType,
            InvalidBitDepthForColorType,
        };
        pub fn unwrap(result: FromBytesResult) UnwrapError!ChunkDataIHDR {
            const width: u31 = switch (result.width) {
                .ok => |width| width,
                .invalid => return error.InvalidWidth,
            };
            const height: u31 = switch (result.height) {
                .ok => |height| height,
                .invalid => return error.InvalidHeight,
            };

            const bit_depth: BitDepth = switch (result.bit_depth) {
                .ok => |bit_depth| bit_depth,
                .invalid => return error.InvalidBitDepth,
                .invalid_for_color => |bit_depth| {
                    std.debug.assert(result.color_type == .ok);
                    std.debug.assert(!result.color_type.ok.bitDepthEnabled(bit_depth));
                    return error.InvalidBitDepthForColorType;
                },
            };
            const color_type: ColorType = switch (result.color_type) {
                .ok => |color_type| color_type,
                .invalid => return error.InvalidColorType,
            };
            std.debug.assert(color_type.bitDepthEnabled(bit_depth));

            return ChunkDataIHDR{
                .width = width,
                .height = height,

                .bit_depth = bit_depth,
                .color_type = color_type,

                .compression_method = result.compression_method,
                .filter_method = result.filter_method,
                .interlace_method = result.interlace_method,
            };
        }
    };
    pub fn fromBytes(bytes: [13]u8) FromBytesResult {
        _ = bytes;
        return std.debug.todo("");
    }
};
