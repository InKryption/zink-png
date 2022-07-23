const std = @import("std");
const chunk = @import("chunk.zig");

const IHDR = @This();
width: u31 align(@alignOf(u32)),
height: u31 align(@alignOf(u32)),

/// A valid bit depth & color type combination.
bd_ct: BitDepthColorType align(@alignOf(u16)),

compression_method: CompressionMethod align(@alignOf(u8)),
filter_method: FilterMethod align(@alignOf(u8)),
interlace_method: InterlaceMethod align(@alignOf(u8)),

pub const BitDepth = enum(u8) {
    bd1 = 1,
    bd2 = 2,
    bd4 = 4,
    bd8 = 8,
    bd16 = 16,

    pub fn colorTypeAllowed(bit_depth: BitDepth, color_type: ColorType) bool {
        return color_type.bitDepthAllowed(bit_depth);
    }
};
pub const ColorType = enum(u8) {
    // used:                | palette(1) | color(2) | alpha(4) | value
    // zig fmt: off
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
                .bd1, .bd2, .bd4, .bd8, .bd16 => true,
            },
            .rgb => switch (bit_depth) {
                .bd1, .bd2, .bd4 => false,
                .bd8, .bd16 => true,
            },
            .palette => switch (bit_depth) {
                .bd1, .bd2, .bd4, .bd8 => true,
                .bd16 => false,
            },
            .grayscale_alpha => switch (bit_depth) {
                .bd1, .bd2, .bd4 => false,
                .bd8, .bd16 => true,
            },
            .rgb_alpha => switch (bit_depth) {
                .bd1, .bd2, .bd4 => false,
                .bd8, .bd16 => true,
            },
        };
    }
};
pub const BitDepthColorType = enum(u16) {
    // -- grayscale --
    bd1_grayscale = enumMemberValue(.bd1, .grayscale),
    bd2_grayscale = enumMemberValue(.bd2, .grayscale),
    bd4_grayscale = enumMemberValue(.bd4, .grayscale),
    bd8_grayscale = enumMemberValue(.bd8, .grayscale),
    bd16_grayscale = enumMemberValue(.bd16, .grayscale),

    // -- rgb --
    bd8_rgb = enumMemberValue(.bd8, .rgb),
    bd16_rgb = enumMemberValue(.bd16, .rgb),

    // -- palette --
    bd1_palette = enumMemberValue(.bd1, .palette),
    bd2_palette = enumMemberValue(.bd2, .palette),
    bd4_palette = enumMemberValue(.bd4, .palette),
    bd8_palette = enumMemberValue(.bd8, .palette),

    // -- grayscale_alpha --
    bd8_grayscale_alpha = enumMemberValue(.bd8, .grayscale_alpha),
    bd16_grayscale_alpha = enumMemberValue(.bd16, .grayscale_alpha),

    // -- rgb_alpha --
    bd8_rgb_alpha = enumMemberValue(.bd8, .rgb_alpha),
    bd16_rgb_alpha = enumMemberValue(.bd16, .rgb_alpha),

    pub fn from(bit_depth: BitDepth, color_type: ColorType) ?BitDepthColorType {
        if (!bit_depth.colorTypeAllowed(color_type)) return null;
        return @intToEnum(BitDepthColorType, enumMemberValue(bit_depth, color_type));
    }

    pub fn bitDepth(bd_ct: BitDepthColorType) BitDepth {
        return bd_ct.values().bd;
    }
    pub fn bitDepthPtr(bd_ct: *BitDepthColorType) *BitDepth {
        return &bd_ct.values().valuesPtr().bd;
    }

    pub fn colorType(bd_ct: BitDepthColorType) ColorType {
        return bd_ct.values().ct;
    }
    pub fn colorTypePtr(bd_ct: *BitDepthColorType) *ColorType {
        return &bd_ct.valuesPtr().ct;
    }

    const Values = packed struct { bd: BitDepth, ct: ColorType };
    fn enumMemberValue(bit_depth: BitDepth, color_type: ColorType) u16 {
        std.debug.assert(bit_depth.colorTypeAllowed(color_type));
        return @bitCast(u16, Values{
            .bd = bit_depth,
            .ct = color_type,
        });
    }
    fn valuesPtr(bd_ct: *BitDepthColorType) *align(@alignOf(BitDepthColorType)) Values {
        return std.mem.bytesAsValue(Values, std.mem.asBytes(bd_ct));
    }
    fn values(bd_ct: BitDepthColorType) Values {
        return @bitCast(Values, bd_ct);
    }

    // Just some assertions that allow for some assumptions.
    comptime {
        const vals = BitDepthColorType.Values{
            .bd = .bd16,
            .ct = .rgb_alpha,
        };
        std.debug.assert(@enumToInt(vals.bd) != @enumToInt(vals.ct));

        const vals_bytes: [2]u8 align(@alignOf(u16)) = @bitCast([2]u8, vals);
        std.debug.assert(vals_bytes[0] == @enumToInt(vals.bd));
        std.debug.assert(vals_bytes[1] == @enumToInt(vals.ct));
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

pub fn toBytes(ihdr: IHDR) [13]u8 {
    var result: [13]u8 = undefined;
    var fbs = std.io.fixedBufferStream(result[0..]);

    fbs.writer().writeIntBig(u32, ihdr.width) catch unreachable;
    fbs.writer().writeIntBig(u32, ihdr.height) catch unreachable;
    fbs.writer().writeIntBig(u8, @enumToInt(ihdr.bd_ct.bitDepth())) catch unreachable;
    fbs.writer().writeIntBig(u8, @enumToInt(ihdr.bd_ct.colorType())) catch unreachable;
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
        /// the original value - (std.math.maxInt(u31) + 1).
        overflow: u31,

        const max = std.math.maxInt(u31);
        /// Returns the real value of the original integer.
        pub fn realValue(ranged: U31RangedResult) u32 {
            return switch (ranged) {
                .ok => |value| value,
                .zero => 0,
                .overflow => |value| @as(u32, value) + (max + 1),
            };
        }
        pub fn fromInt(val: u32) U31RangedResult {
            return switch (val) {
                0 => .zero,
                1...max => .{ .ok = @intCast(u31, val) },
                max + 1...std.math.maxInt(u32) => .{ .overflow = @intCast(u31, val - (max + 1)) },
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
    pub fn unwrap(result: FromBytesResult) UnwrapError!IHDR {
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

        return IHDR{
            .width = width,
            .height = height,

            .bd_ct = BitDepthColorType.from(bit_depth, color_type).?,

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

pub const chunk_type: chunk.Header.Type = .IHDR;
pub fn calculateCrc(ihdr: IHDR) u32 {
    var hasher = std.hash.Crc32.init();
    hasher.hash(&chunk.Header.Type.string(chunk_type));
    hasher.hash(&ihdr.toBytes());
    return hasher.final();
}

test "IHDR" {
    const ihdr = IHDR{
        .width = 100,
        .height = 100,

        .bd_ct = .bd1_grayscale,

        .compression_method = .@"0",
        .filter_method = .@"0",
        .interlace_method = .none,
    };
    try std.testing.expectEqual(ihdr, IHDR.fromBytes(ihdr.toBytes()).unwrap() catch @panic("Something here is definitely wrong"));

    var wrong_bytes align(@alignOf(IHDR)) = ihdr.toBytes();
    std.mem.bytesAsValue(u32, wrong_bytes[0..][0..4]).* = std.math.maxInt(u32);
    std.mem.bytesAsValue(u32, wrong_bytes[4..][0..4]).* = 0;
    std.mem.bytesAsValue(u8, wrong_bytes[8..][0..1]).* = std.math.maxInt(u8);
    std.mem.bytesAsValue(u8, wrong_bytes[9..][0..1]).* = std.math.maxInt(u8);
    std.mem.bytesAsValue(u8, wrong_bytes[10..][0..1]).* = std.math.maxInt(u8);
    std.mem.bytesAsValue(u8, wrong_bytes[11..][0..1]).* = std.math.maxInt(u8);
    std.mem.bytesAsValue(u8, wrong_bytes[12..][0..1]).* = std.math.maxInt(u8);

    try std.testing.expectEqual(IHDR.FromBytesResult{
        .width = .{ .overflow = std.math.maxInt(u32) - (std.math.maxInt(u31) + 1) },
        .height = .zero,

        .bit_depth = .{ .invalid = std.math.maxInt(u8) },
        .color_type = .{ .invalid = std.math.maxInt(u8) },

        .compression_method = @intToEnum(CompressionMethod, std.math.maxInt(u8)),
        .filter_method = @intToEnum(FilterMethod, std.math.maxInt(u8)),
        .interlace_method = @intToEnum(InterlaceMethod, std.math.maxInt(u8)),
    }, IHDR.fromBytes(wrong_bytes));
}
