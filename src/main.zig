const std = @import("std");
const builtin = @import("builtin");
const assert = std.debug.assert;

pub const png_signature = [_]u8{ 137, 80, 78, 71, 13, 10, 26, 10 };

pub const ChunkHeader = struct {
    length: u31,
    type: ChunkType,

    pub const FromBytesError = error{
        InvalidLength,
        InvalidChunkType,
    };
    pub fn fromBytes(bytes: *const [8]u8) FromBytesError!ChunkHeader {
        const length: u31 = try (std.math.cast(u31, std.mem.readIntBig(u32, bytes[0..4])) orelse error.InvalidLength);
        const chunk_type = try (ChunkType.from(bytes[4..8].*) orelse error.InvalidChunkType);
        return ChunkHeader{
            .length = length,
            .type = chunk_type,
        };
    }

    pub fn toBytes(chunk_header: ChunkHeader) [8]u8 {
        var bytes: [8]u8 = undefined;
        std.mem.writeIntBig(u32, bytes[0..4], chunk_header.length);
        bytes[4..8].* = chunk_header.type.string();
        return bytes;
    }
};

test ChunkHeader {
    const header = ChunkHeader{
        .length = 13,
        .type = .IHDR,
    };
    try std.testing.expectEqual(@as(ChunkHeader.FromBytesError!ChunkHeader, header), ChunkHeader.fromBytes(&header.toBytes()));
}

pub const ChunkType = enum(u32) {
    // critical chunks
    IHDR = valueFromString("IHDR".*),
    PLTE = valueFromString("PLTE".*),
    IDAT = valueFromString("IDAT".*),
    IEND = valueFromString("IEND".*),

    // standard ancillary chunks
    tRNS = valueFromString("tRNS".*),
    gAMA = valueFromString("gAMA".*),
    cHRM = valueFromString("cHRM".*),
    sRGB = valueFromString("sRGB".*),
    iCCP = valueFromString("iCCP".*),
    tEXt = valueFromString("tEXt".*),
    zTXt = valueFromString("zTXt".*),
    iTXt = valueFromString("iTXt".*),
    bKGD = valueFromString("bKGD".*),
    pHYs = valueFromString("pHYs".*),
    sBIT = valueFromString("sBIT".*),
    sPLT = valueFromString("sPLT".*),
    hIST = valueFromString("hIST".*),
    tIME = valueFromString("tIME".*),
    _,

    pub inline fn from(bytes: [4]u8) ?ChunkType {
        const ct = @intToEnum(ChunkType, valueFromString(bytes));
        if (!ct.isValidAscii()) return null;
        return ct;
    }

    pub inline fn string(ct: ChunkType) [4]u8 {
        return @bitCast([4]u8, @enumToInt(ct));
    }

    pub inline fn isValidAscii(ct: ChunkType) bool {
        const bytes: [4]u8 = @bitCast([4]u8, @enumToInt(ct));
        inline for (bytes) |byte| {
            if (!std.ascii.isAlphabetic(byte)) return false;
        }
        return true;
    }

    pub inline fn isAncillary(ct: ChunkType) bool {
        return propertyBitMask(ct, 0) != 0;
    }
    pub inline fn isPrivate(ct: ChunkType) bool {
        return propertyBitMask(ct, 1) != 0;
    }
    pub fn isReserved(ct: ChunkType) bool {
        return propertyBitMask(ct, 2) != 0;
    }
    pub inline fn isSafeToCopy(ct: ChunkType) bool {
        return propertyBitMask(ct, 3) != 0;
    }

    inline fn propertyBitMask(ct: ChunkType, index: u2) u8 {
        return string(ct)[index] & 0b00100000;
    }

    inline fn valueFromString(str: [4]u8) u32 {
        return @bitCast(u32, str);
    }
};

test ChunkType {
    try std.testing.expectEqual(@as(?ChunkType, null), ChunkType.from("0ABC".*));
    try std.testing.expect(ChunkType.from("fOOb".*) != null);
    try std.testing.expectEqual(@as(?ChunkType, ChunkType.IHDR), ChunkType.from("IHDR".*));
    inline for (@typeInfo(ChunkType).Enum.fields) |field| {
        try std.testing.expect(field.name.len == 4);
        try std.testing.expect(field.value == @bitCast(u32, field.name[0..4].*));

        const tag: ChunkType = @intToEnum(ChunkType, field.value);
        try std.testing.expect(ChunkType.isValidAscii(tag));
        try std.testing.expect(!ChunkType.isReserved(tag)); // only support non-reserved chunks among the named fields, same as the target PNG version
        try std.testing.expect(!ChunkType.isPrivate(tag)); // only allow for public chunk types among the named fields
    }
}

pub const IHDR = struct {
    width: u31,
    height: u31,
    bit_depth: BitDepth,
    color_type: ColorType,
    compression_method: CompressionMethod,
    filter_method: FilterMethod,
    interlace_method: InterlaceMethod,

    pub const FromBytesError = error{
        InvalidWidth,
        InvalidHeight,
        InvalidBitDepth,
        InvalidColorType,
        InvalidBitDepthForColorType,
        InvalidCompressionMethod,
        InvalidFilterMethod,
    };

    pub fn fromBytes(bytes: *const [13]u8) FromBytesError!IHDR {
        comptime var start = 0;
        defer assert(start == bytes.len);

        const width: u31 = blk: {
            const width: u32 = std.mem.bigToNative(u32, @bitCast(u32, bytes[start..][0..4].*));
            start += 4;
            if (width == 0) return error.InvalidWidth;
            break :blk std.math.cast(u31, width) orelse return error.InvalidWidth;
        };

        const height: u31 = blk: {
            const height: u32 = std.mem.bigToNative(u32, @bitCast(u32, bytes[start..][0..4].*));
            start += 4;
            if (height == 0) return error.InvalidHeight;
            break :blk std.math.cast(u31, height) orelse return error.InvalidHeight;
        };

        const bit_depth: BitDepth = try (util.intToEnum(BitDepth, bytes[start]) orelse error.InvalidBitDepth);
        start += 1;

        const color_type: ColorType = try (util.intToEnum(ColorType, bytes[start]) orelse error.InvalidColorType);
        start += 1;

        if (!color_type.allowsBitDepth(bit_depth)) {
            return error.InvalidBitDepthForColorType;
        }

        const compression_method = try (util.intToEnum(CompressionMethod, bytes[start]) orelse error.InvalidCompressionMethod);
        start += 1;

        const filter_method = try (util.intToEnum(FilterMethod, bytes[start]) orelse error.InvalidFilterMethod);
        start += 1;

        const interlace_method = try (util.intToEnum(InterlaceMethod, bytes[start]) orelse error.InvalidFilterMethod);
        start += 1;

        return IHDR{
            .width = width,
            .height = height,
            .bit_depth = bit_depth,
            .color_type = color_type,
            .compression_method = compression_method,
            .filter_method = filter_method,
            .interlace_method = interlace_method,
        };
    }

    pub fn toBytes(data: IHDR) [13]u8 {
        var result: [13]u8 = undefined;
        result[0..][0..4].* = @bitCast([4]u8, std.mem.nativeToBig(u32, data.width));
        result[4..][0..4].* = @bitCast([4]u8, std.mem.nativeToBig(u32, data.height));
        result[8] = @enumToInt(data.bit_depth);
        result[9] = @enumToInt(data.color_type);
        result[10] = @enumToInt(data.compression_method);
        result[11] = @enumToInt(data.filter_method);
        result[12] = @enumToInt(data.interlace_method);
        return result;
    }

    pub const BitDepth = enum(u8) {
        @"1" = 1,
        @"2" = 2,
        @"4" = 4,
        @"8" = 8,
        @"16" = 16,
    };
    pub const ColorType = enum(u8) {
        grayscale = 0,
        rgb = 2,
        palette = 3,
        grayscale_alpha = 4,
        rgba = 6,

        pub inline fn allowsBitDepth(color_type: ColorType, bit_depth: BitDepth) bool {
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
                .rgba => switch (bit_depth) {
                    .@"1", .@"2", .@"4" => false,
                    .@"8", .@"16" => true,
                },
            };
        }
    };
    pub const CompressionMethod = enum(u8) {
        /// deflate/inflate compression with a sliding window of at most 32768 bytes
        method_0 = 0,
    };
    pub const FilterMethod = enum(u8) {
        /// adaptive filtering with five basic filter types
        method_0 = 0,
    };
    pub const InterlaceMethod = enum(u8) {
        none = 0,
        adam7 = 1,
    };
};
test IHDR {
    const data = IHDR{
        .width = 1,
        .height = 1,
        .bit_depth = .@"1",
        .color_type = .grayscale,
        .compression_method = .method_0,
        .filter_method = .method_0,
        .interlace_method = .none,
    };
    try std.testing.expectEqual(@as(IHDR.FromBytesError!IHDR, data), IHDR.fromBytes(&data.toBytes()));
}

pub const PLTE = struct {
    entries_buf: [256]Entry,
    len_minus_one: LenMinusOne,

    pub const Entry = extern struct { r: u8, g: u8, b: u8 };
    pub const LenMinusOne = std.math.IntFittingRange(1 - 1, 256 - 1);

    pub const FromBytesError = error{
        NoEntries,
        TooManyBytes,
        MissingOneByte,
        MissingTwoBytes,
    };
    pub fn fromBytes(bytes: []const u8) FromBytesError!PLTE {
        if (bytes.len == 0) return error.NoEntries;
        if (bytes.len > 256 * 3) return error.TooManyBytes;
        try switch (bytes.len % 3) {
            0 => {},
            1 => error.MissingTwoBytes,
            2 => error.MissingOneByte,
            else => unreachable,
        };
        const entries_slice = std.mem.bytesAsSlice(Entry, bytes);
        assert(entries_slice.len <= 256);
        var result = PLTE{
            .entries_buf = undefined,
            .len_minus_one = @intCast(LenMinusOne, entries_slice.len - 1),
        };
        std.mem.copy(Entry, result.entries(), entries_slice);
        return result;
    }

    pub fn AsBytes(comptime PointerToPLTE: type) type {
        var pointer_info = @typeInfo(PointerToPLTE).Pointer;
        assert(pointer_info.child == PLTE);
        assert(pointer_info.size == .One);

        pointer_info.child = u8;
        pointer_info.size = .Slice;
        return @Type(.{ .Pointer = pointer_info });
    }
    /// Returns a slice of bytes, backed by the palette entries.
    pub fn asBytes(palette: anytype) AsBytes(@TypeOf(palette)) {
        return std.mem.sliceAsBytes(entries(palette));
    }

    pub inline fn get(palette: *const PLTE, index: LenMinusOne) Entry {
        return palette.entries()[index];
    }

    pub fn Entries(comptime PointerToPLTE: type) type {
        var pointer_info = @typeInfo(PointerToPLTE).Pointer;
        assert(pointer_info.child == PLTE);
        assert(pointer_info.size == .One);

        pointer_info.child = Entry;
        pointer_info.size = .Slice;
        return @Type(.{ .Pointer = pointer_info });
    }
    pub inline fn entries(palette: anytype) Entries(@TypeOf(palette)) {
        return palette.entries_buf[0..len(palette.*)];
    }

    pub const Len = std.math.IntFittingRange(1, 256);
    pub inline fn len(data: PLTE) Len {
        return @as(Len, data.len_minus_one) + 1;
    }
};
test PLTE {
    const data = PLTE{
        .entries_buf = [_]PLTE.Entry{.{ .r = 0, .g = 0, .b = 0 }} ** 256,
        .len_minus_one = 256 - 1,
    };
    try std.testing.expectEqual(@as(PLTE.FromBytesError!PLTE, data), PLTE.fromBytes(data.asBytes()));
}

/// Data for `tRNS` when image is using color type 0
pub const tRNS_0 = enum(u16) {
    _,

    pub inline fn fromBytes(bytes: [2]u8) tRNS_0 {
        const int = std.mem.readIntBig(u16, &bytes);
        return @intToEnum(tRNS_0, int);
    }

    pub inline fn toBytes(gray_level: tRNS_0) [2]u8 {
        var bytes: [2]u8 = undefined;
        std.mem.writeIntBig(u16, &bytes, @enumToInt(gray_level));
        return bytes;
    }

    pub fn WithBitDepth(comptime bit_depth: IHDR.BitDepth) type {
        return switch (bit_depth) {
            .@"1" => u1,
            .@"2" => u2,
            .@"4" => u4,
            .@"8" => u8,
            .@"16" => u16,
        };
    }
    pub inline fn withBitDepth(transparency: tRNS_0, comptime bit_depth: IHDR.BitDepth) WithBitDepth(bit_depth) {
        const int = @enumToInt(transparency);
        const result = @truncate(WithBitDepth(bit_depth), int);
        assert(int == result); // the discarded most significant bits should be 0
        return result;
    }
};

/// Data for `tRNS` when image is using color type 2
pub const tRNS_2 = struct {
    r_level: u16,
    g_level: u16,
    b_level: u16,

    pub inline fn fromBytes(bytes: [6]u8) tRNS_2 {
        return tRNS_2{
            .r_level = std.mem.readIntBig(u16, bytes[0..2]),
            .g_level = std.mem.readIntBig(u16, bytes[2..4]),
            .b_level = std.mem.readIntBig(u16, bytes[4..6]),
        };
    }

    pub inline fn toBytes(transparency: tRNS_2) [6]u8 {
        var bytes: [6]u8 = undefined;
        std.mem.writeIntBig(u16, bytes[0..2], transparency.r_level);
        std.mem.writeIntBig(u16, bytes[2..4], transparency.g_level);
        std.mem.writeIntBig(u16, bytes[4..6], transparency.b_level);
        return bytes;
    }

    pub fn WithBitDepth(comptime bit_depth: IHDR.BitDepth) type {
        const Int = switch (bit_depth) {
            .@"1" => u1,
            .@"2" => u2,
            .@"4" => u4,
            .@"8" => u8,
            .@"16" => u16,
        };
        return struct { r: Int, g: Int, b: Int };
    }
    pub inline fn withBitDepth(transparency: tRNS_2, comptime bit_depth: IHDR.BitDepth) WithBitDepth(bit_depth) {
        const FieldType = std.meta.FieldType;
        const Color = WithBitDepth(bit_depth);
        const result = Color{
            .r = @truncate(FieldType(Color, .r), transparency.r_level),
            .g = @truncate(FieldType(Color, .g), transparency.g_level),
            .b = @truncate(FieldType(Color, .b), transparency.b_level),
        };
        assert(result.r == transparency.r_level); // the discarded most significant bits should be 0
        assert(result.g == transparency.g_level); // the discarded most significant bits should be 0
        assert(result.b == transparency.b_level); // the discarded most significant bits should be 0
        return result;
    }
};

/// Data for `tRNS` when image is using color type 3
pub const tRNS_3 = struct {
    entries_buf: [256]Entry,
    len_minus_one: LenMinusOne,

    pub const Entry = u8;
    pub const LenMinusOne = std.math.IntFittingRange(1 - 1, 256 - 1);

    pub const FromBytesError = error{
        NoEntries,
        TooManyBytes,
    };
    pub fn fromBytes(bytes: []const u8) FromBytesError!tRNS_3 {
        if (bytes.len == 0) return error.NoEntries;
        if (bytes.len > 256) return error.TooManyBytes;
        var result = tRNS_3{
            .entries_buf = undefined,
            .len_minus_one = @intCast(LenMinusOne, bytes.len - 1),
        };
        std.mem.copy(Entry, result.entries(), bytes);
        return result;
    }

    pub fn AsBytes(comptime PointerToPLTE: type) type {
        var pointer_info = @typeInfo(PointerToPLTE).Pointer;
        assert(pointer_info.child == tRNS_3);
        assert(pointer_info.size == .One);

        pointer_info.alignment = 1;
        pointer_info.child = u8;
        pointer_info.size = .Slice;
        return @Type(.{ .Pointer = pointer_info });
    }
    /// Returns a slice of bytes, backed by the transparency entries.
    pub fn asBytes(transparency: anytype) AsBytes(@TypeOf(transparency)) {
        return std.mem.sliceAsBytes(entries(transparency));
    }

    pub inline fn get(transparency: *const tRNS_3, index: LenMinusOne) Entry {
        return transparency.entries()[index];
    }

    pub fn Entries(comptime PointerToPLTE: type) type {
        var pointer_info = @typeInfo(PointerToPLTE).Pointer;
        assert(pointer_info.child == tRNS_3);
        assert(pointer_info.size == .One);

        pointer_info.child = Entry;
        pointer_info.size = .Slice;
        return @Type(.{ .Pointer = pointer_info });
    }
    pub inline fn entries(transparency: anytype) Entries(@TypeOf(transparency)) {
        return transparency.entries_buf[0..len(transparency.*)];
    }

    pub const Len = std.math.IntFittingRange(1, 256);
    pub inline fn len(data: tRNS_3) Len {
        return @as(Len, data.len_minus_one) + 1;
    }
};

pub const TransparentPaletteIterator = struct {
    palette: *const PLTE,
    transparency: *const tRNS_3,
    index: PLTE.LenMinusOne = 0,

    pub const Entry = struct {
        r: u8,
        g: u8,
        b: u8,
        a: u8,
    };

    pub fn next(iter: *TransparentPaletteIterator) ?Entry {
        if (iter.index == iter.palette.len()) return null;
        const rgb = iter.palette.get(iter.index);
        var result = Entry{
            .r = rgb.r,
            .g = rgb.g,
            .b = rgb.b,
            .a = 255,
        };
        if (iter.index < iter.transparency.len()) {
            result.a = iter.transparency.get(iter.index);
        }
        iter.index += 1;
        return result;
    }
};
pub fn transparentPaletteIterator(
    palette: *const PLTE,
    transparency: *const tRNS_3,
) error{TooManyTransparencyEntries}!TransparentPaletteIterator {
    if (palette.len() < transparency.len())
        return error.TooManyTransparencyEntries;
    return TransparentPaletteIterator{
        .palette = palette,
        .transparency = transparency,
    };
}

test transparentPaletteIterator {
    const palette = try PLTE.fromBytes(&[_]u8{
        0x00, 0x00, 0xff,
        0xff, 0xff, 0x00,
        0x32, 0x32, 0x32,
        0xff, 0x00, 0x0,
    });
    const transparency = try tRNS_3.fromBytes(&[_]u8{
        0xef,
        0xab,
        0xfa,
    });
    var iter = try transparentPaletteIterator(&palette, &transparency);
    const Expected = ?TransparentPaletteIterator.Entry;
    try std.testing.expectEqual(@as(Expected, .{ .r = 0x00, .g = 0x00, .b = 0xff, .a = 0xef }), iter.next());
    try std.testing.expectEqual(@as(Expected, .{ .r = 0xff, .g = 0xff, .b = 0x00, .a = 0xab }), iter.next());
    try std.testing.expectEqual(@as(Expected, .{ .r = 0x32, .g = 0x32, .b = 0x32, .a = 0xfa }), iter.next());
    try std.testing.expectEqual(@as(Expected, .{ .r = 0xff, .g = 0x00, .b = 0x00, .a = 0xff }), iter.next());
    try std.testing.expectEqual(@as(Expected, null), iter.next());
    try std.testing.expectEqual(@as(Expected, null), iter.next());
}

/// A typed integer representing the direct value of the `gAMA` chunk.
pub const gAMA = enum(u32) {
    _,

    pub inline fn fromBytes(bytes: [4]u8) gAMA {
        const int = std.mem.readIntBig(u32, &bytes);
        return gAMA.fromInt(int);
    }

    pub inline fn toBytes(gamma: gAMA) [4]u8 {
        var bytes: [4]u8 = undefined;
        const int: u32 = @enumToInt(gamma);
        std.mem.writeIntBig(u32, &bytes, int);
        return bytes;
    }

    pub inline fn fromInt(int_times_100k: u32) gAMA {
        return @intToEnum(gAMA, int_times_100k);
    }

    pub inline fn times_100_000(gamma: gAMA) u32 {
        return @enumToInt(gamma);
    }

    /// The actual gamma value in floating point representation.
    pub inline fn as(gamma: gAMA, comptime T: type) T {
        return @intToFloat(T, gamma.times_100_000()) / 100_000.0;
    }
};

test gAMA {
    const gamma: gAMA = sRGB.compat_gAMA();
    try std.testing.expectEqual(@as(f128, 0.45_455), gamma.as(f128));
    try std.testing.expectEqual(@as(u32, 45_455), gamma.times_100_000());
    try std.testing.expectEqual(gamma, gAMA.fromBytes(gamma.toBytes()));
}

pub const cHRM = struct {
    white_point_x: ValueTimes100K,
    white_point_y: ValueTimes100K,
    red_x: ValueTimes100K,
    red_y: ValueTimes100K,
    green_x: ValueTimes100K,
    green_y: ValueTimes100K,
    blue_x: ValueTimes100K,
    blue_y: ValueTimes100K,

    pub fn fromBytes(bytes: *const [32]u8) cHRM {
        return cHRM{
            // zig fmt: off
            .white_point_x = ValueTimes100K.fromInt(std.mem.readIntBig(u32, bytes[ 0.. 4])),
            .white_point_y = ValueTimes100K.fromInt(std.mem.readIntBig(u32, bytes[ 4.. 8])),
            .red_x =         ValueTimes100K.fromInt(std.mem.readIntBig(u32, bytes[ 8..12])),
            .red_y =         ValueTimes100K.fromInt(std.mem.readIntBig(u32, bytes[12..16])),
            .green_x =       ValueTimes100K.fromInt(std.mem.readIntBig(u32, bytes[16..20])),
            .green_y =       ValueTimes100K.fromInt(std.mem.readIntBig(u32, bytes[20..24])),
            .blue_x =        ValueTimes100K.fromInt(std.mem.readIntBig(u32, bytes[24..28])),
            .blue_y =        ValueTimes100K.fromInt(std.mem.readIntBig(u32, bytes[28..32])),
            // zig fmt: on
        };
    }

    pub fn toBytes(chroma: cHRM) [32]u8 {
        var bytes: [32]u8 = undefined;
        // zig fmt: off
        std.mem.writeIntBig(u32, bytes[ 0.. 4], chroma.white_point_x.times_100_000());
        std.mem.writeIntBig(u32, bytes[ 4.. 8], chroma.white_point_y.times_100_000());
        std.mem.writeIntBig(u32, bytes[ 8..12], chroma.red_x.times_100_000());
        std.mem.writeIntBig(u32, bytes[12..16], chroma.red_y.times_100_000());
        std.mem.writeIntBig(u32, bytes[16..20], chroma.green_x.times_100_000());
        std.mem.writeIntBig(u32, bytes[20..24], chroma.green_y.times_100_000());
        std.mem.writeIntBig(u32, bytes[24..28], chroma.blue_x.times_100_000());
        std.mem.writeIntBig(u32, bytes[28..32], chroma.blue_y.times_100_000());
        // zig fmt: on
        return bytes;
    }

    pub const ValueTimes100K = enum(u32) {
        _,

        pub inline fn fromInt(int_times_100k: u32) ValueTimes100K {
            return @intToEnum(ValueTimes100K, int_times_100k);
        }

        pub inline fn times_100_000(val: ValueTimes100K) u32 {
            return @enumToInt(val);
        }

        pub inline fn as(val: ValueTimes100K, comptime T: type) T {
            return @intToFloat(T, val.times_100_000()) / 100_000.0;
        }
    };
};

test cHRM {
    const chroma: cHRM = sRGB.compat_cHRM();
    try std.testing.expectEqual(chroma, cHRM.fromBytes(&chroma.toBytes()));

    try std.testing.expectEqual(@as(u32, 31_270), chroma.white_point_x.times_100_000());
    try std.testing.expectEqual(@as(f128, 0.3_127), chroma.white_point_x.as(f128));

    try std.testing.expectEqual(@as(u32, 32_900), chroma.white_point_y.times_100_000());
    try std.testing.expectEqual(@as(f128, 0.329), chroma.white_point_y.as(f128));

    try std.testing.expectEqual(@as(u32, 64_000), chroma.red_x.times_100_000());
    try std.testing.expectEqual(@as(f128, 0.64), chroma.red_x.as(f128));

    try std.testing.expectEqual(@as(u32, 33_000), chroma.red_y.times_100_000());
    try std.testing.expectEqual(@as(f128, 0.33), chroma.red_y.as(f128));

    try std.testing.expectEqual(@as(u32, 30_000), chroma.green_x.times_100_000());
    try std.testing.expectEqual(@as(f128, 0.3), chroma.green_x.as(f128));

    try std.testing.expectEqual(@as(u32, 60_000), chroma.green_y.times_100_000());
    try std.testing.expectEqual(@as(f128, 0.6), chroma.green_y.as(f128));

    try std.testing.expectEqual(@as(u32, 15_000), chroma.blue_x.times_100_000());
    try std.testing.expectEqual(@as(f128, 0.15), chroma.blue_x.as(f128));

    try std.testing.expectEqual(@as(u32, 6_000), chroma.blue_y.times_100_000());
    try std.testing.expectEqual(@as(f128, 0.06), chroma.blue_y.as(f128));
}

/// Represents the rendering intent of the sRGB chunk.
pub const sRGB = enum(u8) {
    perceptual = 0,
    relative_colorimetric = 1,
    saturation = 2,
    absolute_colorimetric = 3,

    pub const FromBytesError = error{InvalidRenderingIntent};
    pub inline fn fromBytes(bytes: [1]u8) FromBytesError!sRGB {
        return util.intToEnum(sRGB, bytes[0]) orelse error.InvalidRenderingIntent;
    }

    pub inline fn toBytes(render_intent: sRGB) [1]u8 {
        return @bitCast([1]u8, @enumToInt(render_intent));
    }

    /// The only `gAMA` chunk value allowed for an application writing both
    /// an `sRGB` chunk and a `gAMA` chunk for compatibility.
    pub inline fn compat_gAMA() gAMA {
        comptime return gAMA.fromInt(45455);
    }

    /// The only `cHRM` chunk value allowed for an application writing both
    /// an `sRGB` chunk and a `cHRM` chunk for compatibility.
    pub inline fn compat_cHRM() cHRM {
        comptime return cHRM{
            .white_point_x = cHRM.ValueTimes100K.fromInt(31270),
            .white_point_y = cHRM.ValueTimes100K.fromInt(32900),
            .red_x = cHRM.ValueTimes100K.fromInt(64000),
            .red_y = cHRM.ValueTimes100K.fromInt(33000),
            .green_x = cHRM.ValueTimes100K.fromInt(30000),
            .green_y = cHRM.ValueTimes100K.fromInt(60000),
            .blue_x = cHRM.ValueTimes100K.fromInt(15000),
            .blue_y = cHRM.ValueTimes100K.fromInt(6000),
        };
    }
};

test sRGB {
    const color_space: sRGB = .absolute_colorimetric;
    try std.testing.expectEqual(@as(sRGB.FromBytesError!sRGB, color_space), sRGB.fromBytes(color_space.toBytes()));
}

/// The data at the start of the `iCCP` chunk,
/// preceding the compressed ICC profile.
pub const iCCP_Start = struct {
    name: ProfileName,
    compression_method: CompressionMethod,

    pub const FromBytesError = error{
        /// The stream returned no profile name bytes.
        NoProfileName,
        /// The stream returned a 0 sentinel immediately.
        EmptyProfileName,
        /// The stream ended before returning the 0 sentinel.
        UnfinishedProfileName,
        /// The stream either failed to return a 0 sentinel,
        /// or the encoded profile name is too long.
        ProfileNameTooLong,
        /// The stream ended before returning the compression.
        MissingCompressionMethod,
        /// The stream returned an invalid compression method value.
        InvalidCompressionMethod,
    };
    pub fn fromBytes(reader: anytype) (@TypeOf(reader).Error || FromBytesError)!iCCP_Start {
        const name: ProfileName = blk: {
            var name: ProfileName = .{
                .buffer = undefined,
                .len_minus_one = 1 - 1,
            };
            switch (try reader.readAll(name.slice()[0..1])) {
                0 => return error.NoProfileName,
                1 => if (name.slice()[name.len() - 1] == 0)
                    return error.EmptyProfileName,
                else => unreachable,
            }
            // NOTE: this is `name.buffer.len + 1` because `.len` doesn't include the 0 sentinel of the array.
            while (name.len() < name.buffer.len + 1) : (name.len_minus_one += 1) {
                const byte: u8 = byte: {
                    var byte: [1]u8 = undefined;
                    switch (try reader.readAll(&byte)) {
                        0 => return @as(FromBytesError!iCCP_Start, error.UnfinishedProfileName),
                        1 => {},
                        else => unreachable,
                    }
                    break :byte @bitCast(u8, byte);
                };
                name.buffer[name.len()] = byte;
                if (byte == 0) break;
            } else return @as(FromBytesError!iCCP_Start, error.ProfileNameTooLong);

            assert(name.slice()[name.len()] == 0);
            break :blk name;
        };

        const compression_method: CompressionMethod = blk: {
            const int: u8 = int: {
                var byte: [1]u8 = undefined;
                switch (try reader.readAll(&byte)) {
                    0 => return @as(FromBytesError!iCCP_Start, error.MissingCompressionMethod),
                    1 => {},
                    else => unreachable,
                }
                break :int @bitCast(u8, byte);
            };
            break :blk util.intToEnum(CompressionMethod, int) orelse
                return @as(FromBytesError!iCCP_Start, error.InvalidCompressionMethod);
        };

        return iCCP_Start{
            .name = name,
            .compression_method = compression_method,
        };
    }

    pub const AsBytes = std.BoundedArray(u8, 81);
    pub fn asBytes(icc_profile: *const iCCP_Start) AsBytes {
        var bytes = std.BoundedArray(u8, 81){};
        bytes.appendSliceAssumeCapacity(icc_profile.name.slice());
        bytes.appendAssumeCapacity(0);
        bytes.appendAssumeCapacity(@enumToInt(icc_profile.compression_method));
        return bytes;
    }

    pub const CompressionMethod = enum(u8) {
        /// zlib datastream with deflate compression
        method_0,
    };

    pub const ProfileName = struct {
        buffer: [79:0]u8,
        len_minus_one: LenMinusOne,

        pub const LenMinusOne = std.math.IntFittingRange(1 - 1, 79 - 1);

        pub inline fn fromString(string: []const u8) error{ EmptyString, TooLong }!ProfileName {
            var result = ProfileName{
                .buffer = undefined,
                .len_minus_one = undefined,
            };
            if (string.len == 0) return error.EmptyString;
            if (string.len > result.buffer.len) return error.TooLong;
            result.len_minus_one = @intCast(ProfileName.LenMinusOne, string.len - 1);
            result.buffer[result.len()] = 0;
            std.mem.copy(u8, result.slice(), string);
            return result;
        }

        pub fn Slice(comptime PointerToProfileName: type) type {
            var pointer_info = @typeInfo(PointerToProfileName).Pointer;
            assert(pointer_info.child == ProfileName);
            assert(pointer_info.size == .One);
            assert(pointer_info.sentinel == null);

            pointer_info.alignment = 1;
            pointer_info.child = u8;
            pointer_info.size = .Slice;
            pointer_info.sentinel = &@as(u8, 0);
            return @Type(.{ .Pointer = pointer_info });
        }
        pub inline fn slice(name: anytype) Slice(@TypeOf(name)) {
            return name.buffer[0..len(name.*) :0];
        }

        pub const Len = std.math.IntFittingRange(1, 79);
        pub inline fn len(name: ProfileName) Len {
            return @as(Len, name.len_minus_one) + 1;
        }
    };
};

test iCCP_Start {
    const icc_profile = comptime iCCP_Start{
        .name = iCCP_Start.ProfileName.fromString("foobar") catch unreachable,
        .compression_method = .method_0,
    };

    var fbs = std.io.fixedBufferStream(comptime icc_profile.asBytes().constSlice());
    try std.testing.expectEqual(
        @as(iCCP_Start.FromBytesError!iCCP_Start, icc_profile),
        iCCP_Start.fromBytes(fbs.reader()),
    );
    try std.testing.expectEqualStrings("foobar", icc_profile.name.slice());
}

const util = struct {
    inline fn intToEnum(comptime Enum: type, int: @typeInfo(Enum).Enum.tag_type) ?Enum {
        inline for (@typeInfo(Enum).Enum.fields) |field| {
            if (int == field.value)
                return @field(Enum, field.name);
        }
        return null;
    }

    fn MemoizedErrSet(comptime ErrSet: type) type {
        const set = @typeInfo(ErrSet).ErrorSet orelse return anyerror;
        var errors: [set.len]anyerror = undefined;
        for (errors) |*err, i|
            err.* = @field(ErrSet, set[i].name);
        return MemoizedErrSetImpl(errors);
    }

    fn MemoizedErrSetImpl(
        comptime error_count: comptime_int,
        comptime errors: [error_count]anyerror,
    ) type {
        var set: [errors.len]std.builtin.Type.Error = undefined;
        for (set) |*err, i|
            err.* = .{ .name = @errorName(errors[i]) };
        return @Type(.{ .ErrorSet = &errors });
    }
};
