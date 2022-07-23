const std = @import("std");

const IHDR = @This();
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
    pub fn bitDepthPtr(bd_and_ct: *BitDepthColorType) *BitDepth {
        return &std.mem.bytesAsValue(Value, std.mem.asBytes(bd_and_ct)).bit_depth;
    }

    pub fn colorType(bd_and_ct: BitDepthColorType) ColorType {
        return @bitCast(BitDepthColorType.Value, @enumToInt(bd_and_ct)).color_type;
    }
    pub fn colorTypePtr(bd_and_ct: *BitDepthColorType) *BitDepth {
        return &std.mem.bytesAsValue(Value, std.mem.asBytes(bd_and_ct)).color_type;
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

pub fn toBytes(ihdr: IHDR) [13]u8 {
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

            .bit_depth_color_type = BitDepthColorType.from(bit_depth, color_type).?,

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

pub fn bitDepth(ihdr: IHDR) BitDepth {
    return ihdr.bit_depth_color_type.bitDepth();
}
pub fn bitDepthPtr(ihdr: *IHDR) *BitDepth {
    return ihdr.bit_depth_color_type.bitDepthPtr();
}
pub fn colorType(ihdr: IHDR) ColorType {
    return ihdr.bit_depth_color_type.colorType();
}
pub fn colorTypePtr(ihdr: *IHDR) *ColorType {
    ihdr.bit_depth_color_type.colorTypePtr();
}

test "IHDR" {
    const ihdr = IHDR{
        .width = 100,
        .height = 100,

        .bit_depth_color_type = .@"1_grayscale",

        .compression_method = .@"0",
        .filter_method = .@"0",
        .interlace_method = .none,
    };
    try std.testing.expectEqual(ihdr, IHDR.fromBytes(ihdr.toBytes()).unwrap() catch @panic("Something here is definitely wrong"));

    var wrong_bytes = ihdr.toBytes();
    const width: *align(@alignOf(u8)) u32 = std.mem.bytesAsValue(u32, wrong_bytes[0..][0..4]);
    const height: *align(@alignOf(u8)) u32 = std.mem.bytesAsValue(u32, wrong_bytes[4..][0..4]);
    const bit_depth: *align(@alignOf(u8)) u8 = std.mem.bytesAsValue(u8, wrong_bytes[8..][0..1]);
    const color_type: *align(@alignOf(u8)) u8 = std.mem.bytesAsValue(u8, wrong_bytes[9..][0..1]);

    width.* = 0;
    try std.testing.expectError(error.InvalidWidth, IHDR.fromBytes(wrong_bytes).unwrap());

    width.* = std.mem.nativeToBig(u32, ihdr.width);
    height.* = 0;
    try std.testing.expectError(error.InvalidHeight, IHDR.fromBytes(wrong_bytes).unwrap());

    height.* = std.mem.nativeToBig(u32, ihdr.height);
    bit_depth.* = std.mem.nativeToBig(u8, std.math.maxInt(u8));
    try std.testing.expectError(error.InvalidBitDepth, IHDR.fromBytes(wrong_bytes).unwrap());

    bit_depth.* = std.mem.nativeToBig(u8, @enumToInt(ihdr.bitDepth()));
    color_type.* = std.mem.nativeToBig(u8, std.math.maxInt(u8));
    try std.testing.expectError(error.InvalidColorType, IHDR.fromBytes(wrong_bytes).unwrap());

    try std.testing.expectEqual(
        @as(?IHDR.BitDepthColorType, null),
        IHDR.BitDepthColorType.from(ihdr.bitDepth(), .rgb),
    );

    color_type.* = std.mem.nativeToBig(u8, @enumToInt(IHDR.ColorType.rgb));
    try std.testing.expectError(error.InvalidColorTypeForBitDepth, IHDR.fromBytes(wrong_bytes).unwrap());

    color_type.* = std.mem.nativeToBig(u8, @enumToInt(ihdr.colorType()));
    try std.testing.expectEqual(ihdr, IHDR.fromBytes(wrong_bytes).unwrap() catch @panic("Welp, something wrong here."));
}
