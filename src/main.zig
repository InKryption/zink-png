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

    pub fn from(bytes: [4]u8) ChunkType {
        return @intToEnum(ChunkType, (std.mem.readIntBig(u32, bytes[0..])));
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
        var fbs = std.io.fixedBufferStream(result[0..]);
        fbs.writer().writeIntBig(u32, header.length) catch unreachable;
        fbs.writer().writeAll(header.type.string()[0..]) catch unreachable;
        return result;
    }

    pub const FromBytesResult = struct {
        length: Length,
        type: Type,

        pub const Length = union(enum) { ok: u31, invalid: u32 };
        /// `ok.isValidAscii()` => `true`
        /// `invalid.isValidAscii()` => `false`
        pub const Type = union(enum) { ok: ChunkType, invalid: ChunkType };

        pub const UnwrapError = error{
            InvalidLength,
            InvalidType,
        };
        pub fn unwrap(result: FromBytesResult) UnwrapError!ChunkHeader {
            return ChunkHeader{
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
                const ch_type = ChunkType.from(fbs.reader().readBytesNoEof(4) catch unreachable);
                break :ch_type if (ch_type.isValidAscii())
                    FromBytesResult.Type{ .ok = ch_type }
                else
                    FromBytesResult.Type{ .invalid = ch_type };
            },
        };
    }
};

test "ChunkType & ChunkHeader" {
    try std.testing.expectEqualStrings("WOAH", &ChunkType.from("WOAH".*).string());
    for (@as([4096]void, undefined)) |_, sample_valid_ch_len| {
        inline for (comptime std.enums.values(ChunkType)) |sample_valid_ch_ty| {
            const valid_header = ChunkHeader{
                .length = @intCast(u31, sample_valid_ch_len),
                .type = sample_valid_ch_ty,
            };
            try std.testing.expectEqual(valid_header, ChunkHeader.fromBytes(valid_header.toBytes()).unwrap() catch @panic("That shouldn't have happened."));

            const invalid_len: u32 = std.math.maxInt(u31) + 1 + @intCast(u32, sample_valid_ch_len);
            try std.testing.expectEqual(ChunkHeader.FromBytesResult{
                .length = .{ .invalid = invalid_len },
                .type = .{ .ok = sample_valid_ch_ty },
            }, ChunkHeader.fromBytes(blk: {
                var bytes_with_invalid_len: [8]u8 = undefined;
                bytes_with_invalid_len[0..4].* = std.mem.toBytes(std.mem.nativeToBig(u32, invalid_len));
                bytes_with_invalid_len[4..8].* = sample_valid_ch_ty.string();
                break :blk bytes_with_invalid_len;
            }));
        }
        try std.testing.expectEqual(
            ChunkHeader.FromBytesResult{
                .length = .{ .ok = @intCast(u31, sample_valid_ch_len) },
                .type = .{ .invalid = ChunkType.from("0000".*) },
            },
            ChunkHeader.fromBytes(blk: {
                var bytes_with_invalid_type: [8]u8 = undefined;
                bytes_with_invalid_type[0..4].* = std.mem.toBytes(std.mem.nativeToBig(u32, @intCast(u31, sample_valid_ch_len)));
                bytes_with_invalid_type[4..8].* = ChunkType.from("0000".*).string();
                break :blk bytes_with_invalid_type;
            }),
        );
    }
}

pub const ChunkDataIHDR = struct {
    width: u31 align(@alignOf(u32)),
    height: u31 align(@alignOf(u32)),

    bit_depth_color_type: BitDepthColorType align(@alignOf(u16)),

    compression_method: CompressionMethod align(@alignOf(u8)),
    filter_method: FilterMethod align(@alignOf(u8)),
    interlace_method: InterlaceMethod align(@alignOf(u8)),

    pub const BitDepth = enum(u8) {
        @"1" = 1,
        @"2" = 2,
        @"4" = 4,
        @"8" = 8,
        @"16" = 16,

        pub fn colorTypeAllowed(bit_depth: BitDepth, color_type: ColorType) bool {
            return color_type.bitDepthAllowed(bit_depth);
        }
    };
    pub const ColorType = enum(u8) {
        // zig fmt: off
        // used:                | palette(1) | color(2) | alpha(4) | value
        grayscale               =         0  +       0  +       0,      // 0
        rgb                     =         0  +       2  +       0,      // 2
        palette                 =         1  +       2  +       0,      // 3
        grayscale_alpha         =         0  +       0  +       4,      // 4
        rgb_alpha               =         0  +       2  +       4,      // 6
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

        pub fn bitDepthAllowed(color_type: ColorType, bit_depth: BitDepth) bool {
            return switch (color_type) {
                .grayscale => switch (bit_depth) {
                    .@"1", .@"2", .@"4", .@"8", .@"16" => true,
                },
                .rgb => switch (bit_depth) {
                    .@"1", .@"2", .@"4" => false,
                    .@"8", .@"16" => true,
                },
                .palette => switch (bit_depth) {
                    .@"1", .@"2", .@"4", .@"8" => true,
                    .@"16" => false,
                },
                .grayscale_alpha => switch (bit_depth) {
                    .@"1", .@"2", .@"4" => false,
                    .@"8", .@"16" => true,
                },
                .rgb_alpha => switch (bit_depth) {
                    .@"1", .@"2", .@"4" => false,
                    .@"8", .@"16" => true,
                },
            };
        }
    };
    pub const BitDepthColorType = enum(u16) {
        // -- grayscale --
        @"1_grayscale" = enumValue(.@"1", .grayscale),
        @"2_grayscale" = enumValue(.@"2", .grayscale),
        @"4_grayscale" = enumValue(.@"4", .grayscale),
        @"8_grayscale" = enumValue(.@"8", .grayscale),
        @"16_grayscale" = enumValue(.@"16", .grayscale),

        // -- rgb --
        @"8_rgb" = enumValue(.@"8", .rgb),
        @"16_rgb" = enumValue(.@"16", .rgb),

        // -- palette --
        @"1_palette" = enumValue(.@"1", .palette),
        @"2_palette" = enumValue(.@"2", .palette),
        @"4_palette" = enumValue(.@"4", .palette),
        @"8_palette" = enumValue(.@"8", .palette),

        // -- grayscale_alpha --
        @"8_grayscale_alpha" = enumValue(.@"8", .grayscale_alpha),
        @"16_grayscale_alpha" = enumValue(.@"16", .grayscale_alpha),

        // -- rgb_alpha --
        @"8_rgb_alpha" = enumValue(.@"8", .rgb_alpha),
        @"16_rgb_alpha" = enumValue(.@"16", .rgb_alpha),

        pub fn from(bit_depth: BitDepth, color_type: ColorType) ?BitDepthColorType {
            if (!bit_depth.colorTypeAllowed(color_type)) return null;
            return @intToEnum(BitDepthColorType, enumValue(bit_depth, color_type));
        }

        pub fn bitDepth(bd_and_ct: BitDepthColorType) BitDepth {
            return @bitCast(BitDepthColorType.Value, @enumToInt(bd_and_ct)).bit_depth;
        }

        pub fn colorType(bd_and_ct: BitDepthColorType) ColorType {
            return @bitCast(BitDepthColorType.Value, @enumToInt(bd_and_ct)).color_type;
        }

        const Value = extern struct {
            bit_depth: BitDepth align(@alignOf(u8)),
            color_type: ColorType align(@alignOf(u8)),

            fn toInt(color_type_and_bit_depth_value_bits: Value) u16 {
                return @bitCast(u16, color_type_and_bit_depth_value_bits);
            }
        };
        fn enumValue(bit_depth: BitDepth, color_type: ColorType) u16 {
            std.debug.assert(bit_depth.colorTypeAllowed(color_type));
            const value = Value{
                .bit_depth = bit_depth,
                .color_type = color_type,
            };
            return value.toInt();
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
        var result: [13]u8 = undefined;
        var fbs = std.io.fixedBufferStream(result[0..]);

        fbs.writer().writeIntBig(u32, ihdr.width) catch unreachable;
        fbs.writer().writeIntBig(u32, ihdr.height) catch unreachable;
        fbs.writer().writeIntBig(u8, @enumToInt(ihdr.bitDepth())) catch unreachable;
        fbs.writer().writeIntBig(u8, @enumToInt(ihdr.colorType())) catch unreachable;
        fbs.writer().writeIntBig(u8, @enumToInt(ihdr.compression_method)) catch unreachable;
        fbs.writer().writeIntBig(u8, @enumToInt(ihdr.filter_method)) catch unreachable;
        fbs.writer().writeIntBig(u8, @enumToInt(ihdr.interlace_method)) catch unreachable;
        std.debug.assert(fbs.pos == result.len);

        return result;
    }

    pub const FromBytesResult = struct {
        width: WidthResult,
        height: HeightResult,

        bit_depth: BitDepthResult,
        color_type: ColorTypeResult,

        compression_method: CompressionMethod,
        filter_method: FilterMethod,
        interlace_method: InterlaceMethod,

        const U31RangedResult = union(enum) {
            /// valid value. Returned as is.
            ok: u31,
            /// value was 0. Returned value is implicitly 0.
            zero,
            /// value was > std.math.maxInt(u31). Returned value is
            /// the original value minus std.math.maxInt(u31).
            overflow: u31,
            /// Returns the real value of the original integer.
            pub fn realValue(ranged: U31RangedResult) u32 {
                return switch (ranged) {
                    .ok => |value| value,
                    .zero => 0,
                    .overflow => |value| @as(u32, value) + std.math.maxInt(u31),
                };
            }
            pub fn fromInt(val: u32) U31RangedResult {
                const max = std.math.maxInt(u31);
                return switch (val) {
                    0 => .zero,
                    1...max => .{ .ok = @intCast(u31, val) },
                    max + 1...std.math.maxInt(u32) => .{ .overflow = @intCast(u31, val - max) },
                };
            }
        };
        pub const WidthResult = U31RangedResult;
        pub const HeightResult = U31RangedResult;

        pub const BitDepthResult = union(enum) { ok: BitDepth, invalid: u8 };
        pub const ColorTypeResult = union(enum) { ok: ColorType, invalid_for_bit_depth: ColorType, invalid: u8 };

        pub const UnwrapError = error{
            InvalidWidth,
            InvalidHeight,

            InvalidBitDepth,
            InvalidColorType,
            InvalidColorTypeForBitDepth,
        };
        pub fn unwrap(result: FromBytesResult) UnwrapError!ChunkDataIHDR {
            const width: u31 = switch (result.width) {
                .ok => |width| width,
                .zero, .overflow => return error.InvalidWidth,
            };
            std.debug.assert(width > 0);

            const height: u31 = switch (result.height) {
                .ok => |height| height,
                .zero, .overflow => return error.InvalidHeight,
            };
            std.debug.assert(height > 0);

            const bit_depth: BitDepth = switch (result.bit_depth) {
                .ok => |bit_depth| bit_depth,
                .invalid => return error.InvalidBitDepth,
            };
            const color_type: ColorType = switch (result.color_type) {
                .ok => |color_type| color_type,
                .invalid_for_bit_depth => |invalid_color| {
                    std.debug.assert(!invalid_color.bitDepthAllowed(bit_depth));
                    return error.InvalidColorTypeForBitDepth;
                },
                .invalid => return error.InvalidColorType,
            };
            std.debug.assert(color_type.bitDepthAllowed(bit_depth));

            const compression_method: CompressionMethod = result.compression_method;
            const filter_method: FilterMethod = result.filter_method;
            const interlace_method: InterlaceMethod = result.interlace_method;

            return ChunkDataIHDR{
                .width = width,
                .height = height,

                .bit_depth_color_type = BitDepthColorType.from(bit_depth, color_type),

                .compression_method = compression_method,
                .filter_method = filter_method,
                .interlace_method = interlace_method,
            };
        }
    };
    pub fn fromBytes(bytes: [13]u8) FromBytesResult {
        var fbs = std.io.fixedBufferStream(bytes[0..]);
        defer std.debug.assert(fbs.pos == bytes.len);

        const width = FromBytesResult.WidthResult.fromInt(fbs.reader().readIntBig(u32) catch unreachable);
        const height = FromBytesResult.HeightResult.fromInt(fbs.reader().readIntBig(u32) catch unreachable);

        const bit_depth: FromBytesResult.BitDepthResult = bit_depth: {
            const bit_depth_raw: u8 = fbs.reader().readIntBig(u8) catch unreachable;
            if (std.meta.intToEnum(BitDepth, bit_depth_raw)) |bit_depth| {
                break :bit_depth .{ .ok = bit_depth };
            } else |err| switch (err) {
                error.InvalidEnumTag => break :bit_depth .{ .invalid = bit_depth_raw },
            }
        };
        const color_type: FromBytesResult.ColorTypeResult = color_type: {
            const color_type_raw = fbs.reader().readIntBig(u8) catch unreachable;
            const color_type = std.meta.intToEnum(ColorType, color_type_raw) catch |err| switch (err) {
                error.InvalidEnumTag => break :color_type .{ .invalid = color_type_raw },
            };
            switch (bit_depth) {
                .ok => |bd| if (!bd.colorTypeAllowed(color_type)) {
                    break :color_type .{ .invalid_for_bit_depth = color_type };
                },
                .invalid => {},
            }
            break :color_type .{ .ok = color_type };
        };

        const compression_method = @intToEnum(CompressionMethod, fbs.reader().readIntBig(u8) catch unreachable);
        const filter_method = @intToEnum(FilterMethod, fbs.reader().readIntBig(u8) catch unreachable);
        const interlace_method = @intToEnum(InterlaceMethod, fbs.reader().readIntBig(u8) catch unreachable);

        return FromBytesResult{
            .width = width,
            .height = height,

            .bit_depth = bit_depth,
            .color_type = color_type,

            .compression_method = compression_method,
            .filter_method = filter_method,
            .interlace_method = interlace_method,
        };
    }

    pub fn bitDepth(ihdr: ChunkDataIHDR) BitDepth {
        return ihdr.bit_depth_color_type.bitDepth();
    }
    pub fn colorType(ihdr: ChunkDataIHDR) ColorType {
        return ihdr.bit_depth_color_type.colorType();
    }
};

test "ChunkDataIHDR" {
    var ihdr = ChunkDataIHDR{
        .width = 100,
        .height = 100,

        .bit_depth_color_type = .@"1_grayscale",

        .compression_method = .@"0",
        .filter_method = .@"0",
        .interlace_method = .none,
    };
    try std.testing.expectEqual(ihdr, ChunkDataIHDR.fromBytes(ihdr.toBytes()).unwrap() catch @panic("Something here is definitely wrong"));
}

pub const ChunkDataPLTE = struct {
    entries_buf: [256]Entry,
    count_minus_one: u8,

    pub const Entry = extern struct { r: u8, g: u8, b: u8 };
    pub fn count(plte: ChunkDataPLTE) std.math.IntFittingRange(1, 256) {
        return @as(std.math.IntFittingRange(1, 256), plte.count_minus_one) + 1;
    }
    pub fn entries(plte: *const ChunkDataPLTE) []const Entry {
        return plte.entries_buf[0..plte.count()];
    }

    pub const RawBytes = std.BoundedArray(u8, 256 * 3);

    pub fn toBytes(plte: ChunkDataPLTE) RawBytes {
        var result = RawBytes{};
        for (plte.entries()) |entry| {
            result.appendSlice(&[_]u8{ entry.r, entry.g, entry.b }) catch unreachable;
        }
        return result;
    }

    pub const FromBytesResult = union(enum) {
        ok: ChunkDataPLTE,
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
    };
    /// Asserts `bytes.len > 0`.
    pub fn fromBytes(bytes: RawBytes) FromBytesResult {
        std.debug.assert(bytes.len > 0);

        var result_entries = std.BoundedArray(Entry, 256){};
        const used_bytes_len = bytes.len - (bytes.len % 3);

        comptime std.debug.assert(@typeInfo(Entry).Struct.layout == .Extern);
        comptime std.debug.assert(std.mem.eql(u8, @typeInfo(Entry).Struct.fields[0].name, "r"));
        comptime std.debug.assert(std.mem.eql(u8, @typeInfo(Entry).Struct.fields[1].name, "g"));
        comptime std.debug.assert(std.mem.eql(u8, @typeInfo(Entry).Struct.fields[2].name, "b"));
        result_entries.appendSliceAssumeCapacity(std.mem.bytesAsSlice(Entry, bytes.constSlice()[0..used_bytes_len]));

        return switch (bytes.len - used_bytes_len) {
            0 => FromBytesResult{ .ok = ChunkDataPLTE{
                .entries_buf = result_entries.buffer,
                .count_minus_one = @intCast(u8, result_entries.len - 1),
            } },
            1 => FromBytesResult{ .missing_one_byte = FromBytesResult.MissingOneByte{
                .full = std.BoundedArray(Entry, 255).fromSlice(entries.constSlice()) catch unreachable,
                .trailing = .{
                    .r = bytes.get(bytes.len - 2),
                    .g = bytes.get(bytes.len - 1),
                },
            } },
            2 => FromBytesResult{ .missing_two_bytes = FromBytesResult.MissingTwoBytes{
                .full = std.BoundedArray(Entry, 255).fromSlice(entries.constSlice()) catch unreachable,
                .trailing = .{
                    .r = bytes.get(bytes.len - 1),
                },
            } },
            else => unreachable,
        };
    }
};

test "ChunkDataPLTE" {
    _ = ChunkDataPLTE;
    std.log.warn("\nTODO: test stuff here.", .{});
    return error.SkipZigTest;
}

pub const ChunkDataIDAT = struct {
    bytes: []const u8,
};

pub const ChunkDataIEND = struct {};
