const std = @import("std");
const builtin = @import("builtin");
const assert = std.debug.assert;

pub const png_signature = [_]u8{ 137, 80, 78, 71, 13, 10, 26, 10 };

pub const ChunkHeader = struct {
    length: u31,
    type: ChunkType,

    pub const ValidateError = error{
        InvalidLength,
        InvalidType,
    };
    pub fn fromBytes(bytes: *const [8]u8) Raw {
        var result: Raw = @bitCast(Raw, bytes.*);
        result.length = std.mem.bigToNative(u32, result.length);
        return result;
    }

    pub fn toRaw(chunk_header: ChunkHeader) Raw {
        return Raw{
            .length = chunk_header.length,
            .type = chunk_header.type.string(),
        };
    }

    pub const Raw = extern struct {
        length: u32,
        type: [4]u8,

        pub fn validate(raw: Raw) ValidateError!ChunkHeader {
            return ChunkHeader{
                .length = try (std.math.cast(u31, raw.length) orelse error.InvalidLength),
                .type = try (ChunkType.from(raw.type) orelse error.InvalidType),
            };
        }

        pub fn toBytes(raw: Raw) [8]u8 {
            var copy = raw;
            copy.length = std.mem.nativeToBig(u32, copy.length);
            return @bitCast([8]u8, copy);
        }
    };
};

test ChunkHeader {
    const header = ChunkHeader{
        .length = 13,
        .type = .IHDR,
    };
    try std.testing.expectEqual(@as(ChunkHeader.ValidateError!ChunkHeader, header), ChunkHeader.fromBytes(&header.toRaw().toBytes()).validate());
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

pub fn ChunkIterator(
    comptime ReaderType: type,
    /// Enables calculating the CRC code.
    comptime calculate_crc: bool,
) type {
    return struct {
        const Self = @This();
        inner: Inner,
        remaining_bytes: u31 = undefined,
        crc_hasher: if (calculate_crc) std.hash.Crc32 else void = if (calculate_crc) std.hash.Crc32.init() else {},
        expected_crc: ?u32 = null,
        state: State = .start,

        const State = enum {
            start,
            end,
            expecting_header,
            exhausting_data,
            got_expected_crc,
        };

        pub const SignatureCheck = ChunkIteratorSignatureCheck;
        /// The result describes whether the PNG
        /// signature was attained correctly from
        /// the stream.
        pub fn checkSignature(iter: *Self) Inner.Error!SignatureCheck {
            var out: SignatureCheck = undefined;
            try @call(.always_inline, checkSignatureAdvanced, .{ iter, &out });
            return out;
        }
        /// The result is written to `out`, whether
        /// returning normally or with an error. It
        /// describes whether the PNG signature was
        /// attained correctly from the stream.
        pub fn checkSignatureAdvanced(iter: *Self, out: *SignatureCheck) Inner.Error!void {
            switch (iter.state) {
                .start => {},
                .end,
                .expecting_header,
                .exhausting_data,
                .got_expected_crc,
                => unreachable,
            }

            out.* = undefined;
            defer switch (out.*) {
                else => {},
            };

            var bytes: std.BoundedArray(u8, png_signature.len) = .{};
            iter.inner.readIntoBoundedBytes(png_signature.len, &bytes) catch |err| {
                assert(bytes.len < png_signature.len);
                out.* = .{ .incomplete_signature = std.BoundedArray(u8, png_signature.len).fromSlice(bytes.constSlice()) catch unreachable };
                iter.state = .end;
                return err;
            };
            if (bytes.len < png_signature.len) {
                out.* = .{ .incomplete_signature = std.BoundedArray(u8, png_signature.len - 1).fromSlice(bytes.constSlice()) catch unreachable };
                iter.state = .end;
                return;
            }
            if (!std.mem.eql(u8, bytes.constSlice(), &png_signature)) {
                out.* = .{ .bad_signature = bytes.constSlice()[0..png_signature.len].* };
                iter.state = .end;
                return;
            }

            out.* = .ok;
            iter.state = .expecting_header;
            return;
        }

        pub const NextResult = ChunkIteratorNextResult;
        /// The result describes whether the chunk
        /// header was attained from the stream and
        /// decoded correctly.
        pub fn next(iter: *Self) Inner.Error!NextResult {
            var out: NextResult = undefined;
            try @call(.always_inline, nextAdvanced, .{ iter, &out });
            return out;
        }
        /// The result is written to `out`, whether
        /// returning normally or with an error. It
        /// describes whether the chunk header was
        /// attained from the stream and decoded correctly.
        pub fn nextAdvanced(iter: *Self, out: *NextResult) Inner.Error!void {
            switch (iter.state) {
                .start => unreachable, // meant to call `checkSignature*` first.s
                .end => {
                    out.* = .none;
                    return;
                },
                .expecting_header => {
                    assert(iter.expected_crc == null);
                },
                .exhausting_data => unreachable, // supposed to call meant to call `fetchExpectedCrc` first.
                .got_expected_crc => {
                    assert(iter.remaining_bytes == 0);
                    assert(iter.expected_crc != null);
                    iter.remaining_bytes = undefined;
                    iter.expected_crc = null;
                },
            }

            var bytes: std.BoundedArray(u8, 8) = .{};
            iter.inner.readIntoBoundedBytes(8, &bytes) catch |err| {
                if (bytes.len != 0) {
                    out.* = .{ .incomplete_chunk_header = std.BoundedArray(u8, 8 - 1).fromSlice(bytes.constSlice()) catch unreachable };
                } else {
                    out.* = .none;
                }
                iter.state = .end;
                return err;
            };
            if (bytes.len < 8) {
                if (bytes.len != 0) {
                    out.* = .{ .incomplete_chunk_header = std.BoundedArray(u8, 8 - 1).fromSlice(bytes.constSlice()) catch unreachable };
                } else {
                    out.* = .none;
                }
                iter.state = .end;
                return;
            }
            const raw_reader = ChunkHeader.fromBytes(bytes.constSlice()[0..8]);
            const header = raw_reader.validate() catch |err| {
                out.* = switch (err) {
                    error.InvalidLength,
                    error.InvalidType,
                    => .{ .invalid_chunk_header = raw_reader },
                };
                iter.state = .end;
                return;
            };
            out.* = .{ .ok = header };
            iter.state = .exhausting_data;
            iter.remaining_bytes = header.length;
            iter.crc_hasher = std.hash.Crc32.init();
            iter.crc_hasher.update(&header.type.string());
        }

        /// Returns true on successfully reading the CRC code,
        /// false if the stream ended before returning it fully.
        /// Returning a read error implies the same as returning false.
        /// Also see `getExpectedCrc` and `getActualCrc`.
        pub fn fetchExpectedCrc(iter: *Self) Inner.Error!bool {
            switch (iter.state) {
                .start,
                .end,
                .expecting_header,
                => unreachable, //  not valid to call yet

                .got_expected_crc => unreachable, // shouldn't call twice
                .exhausting_data => {
                    assert(iter.remaining_bytes == 0); // caller is meant to exhaust the stream returned by `.dataReader()`.
                    assert(iter.expected_crc == null);
                },
            }

            var crc_code_bytes: std.BoundedArray(u8, 4) = .{};
            iter.inner.readIntoBoundedBytes(4, &crc_code_bytes) catch |err| {
                iter.state = .end;
                return err;
            };
            if (crc_code_bytes.len < 4) {
                iter.state = .end;
                return false;
            }
            const crc_code = std.mem.readIntBig(u32, crc_code_bytes.constSlice()[0..4]);
            iter.expected_crc = crc_code;
            iter.state = .got_expected_crc;
            return true;
        }

        /// Get the CRC code found after the chunk data.
        pub inline fn getExpectedCrc(iter: *const Self) u32 {
            assert(iter.remaining_bytes == 0); //
            // if null, the caller probably forgot to call `fetchExpectedCrc`
            return iter.expected_crc.?;
        }

        /// Get the CRC code calculated based off of the
        /// data read from the stream.
        pub inline fn getActualCrc(iter: *const Self) u32 {
            comptime assert(calculate_crc); // must enable CRC calculation in order to use this
            assert(iter.remaining_bytes == 0); // must exhaust the data reader stream first
            var hasher_copy: std.hash.Crc32 = iter.crc_hasher;
            return hasher_copy.final();
        }

        pub const Inner = ReaderType;
        pub const DataReader = std.io.Reader(*Self, Inner.Error, read);
        /// The caller should read all the bytes from this
        /// stream before calling `fetchExpectedCrc`, which
        /// must be called before a subsequent call to `next*`.
        pub inline fn dataReader(iter: *Self) DataReader {
            return .{ .context = iter };
        }

        fn read(iter: *Self, buffer: []u8) Inner.Error!usize {
            if (buffer.len == 0) return 0;
            if (iter.remaining_bytes == 0) return 0;
            const amt = @intCast(u31, @min(buffer.len, iter.remaining_bytes));
            const bytes_read = @intCast(u31, try iter.inner.read(buffer[0..amt]));
            if (calculate_crc) {
                iter.crc_hasher.update(buffer[0..bytes_read]);
            }
            iter.remaining_bytes -= bytes_read;
            return bytes_read;
        }
    };
}
pub inline fn chunkIterator(
    reader: anytype,
    /// Enables calculating the CRC code.
    comptime calculate_crc: bool,
) ChunkIterator(@TypeOf(reader), calculate_crc) {
    return .{ .inner = reader };
}
pub const ChunkIteratorSignatureCheck = union(enum) {
    ok,
    bad_signature: [png_signature.len]u8,
    incomplete_signature: std.BoundedArray(u8, png_signature.len - 1),

    pub inline fn unwrap(result: ChunkIteratorSignatureCheck) error{ BadPngSignature, IncompletePngSignature }!void {
        return switch (result) {
            .ok => {},
            .bad_signature => error.BadPngSignature,
            .incomplete_signature => error.IncompletePngSignature,
        };
    }
};
pub const ChunkIteratorNextResult = union(enum) {
    none,
    ok: ChunkHeader,
    incomplete_chunk_header: std.BoundedArray(u8, 8 - 1),
    invalid_chunk_header: ChunkHeader.Raw,

    pub inline fn unwrap(result: ChunkIteratorNextResult) error{ IncompleteBytes, InvalidChunkHeader }!?ChunkHeader {
        return switch (result) {
            .none => null,
            .ok => |header| header,
            .incomplete_chunk_header => error.IncompleteBytes,
            .invalid_chunk_header => error.InvalidChunkHeader,
        };
    }
};

pub const IHDR = struct {
    pub const Data = struct {
        width: u31,
        height: u31,
        bit_depth: BitDepth,
        color_type: ColorType,
        compression_method: CompressionMethod,
        filter_method: FilterMethod,
        interlace_method: InterlaceMethod,

        pub fn toBytes(data: IHDR.Data) [13]u8 {
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
    };

    pub const FromBytesError = error{
        InvalidWidth,
        InvalidHeight,
        InvalidBitDepth,
        InvalidColorType,
        InvalidBitDepthForColorType,
        InvalidCompressionMethod,
        InvalidFilterMethod,
    };

    pub fn fromBytes(bytes: *const [13]u8) FromBytesError!IHDR.Data {
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

        return IHDR.Data{
            .width = width,
            .height = height,
            .bit_depth = bit_depth,
            .color_type = color_type,
            .compression_method = compression_method,
            .filter_method = filter_method,
            .interlace_method = interlace_method,
        };
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
    const data = IHDR.Data{
        .width = 1,
        .height = 1,
        .bit_depth = .@"1",
        .color_type = .grayscale,
        .compression_method = .method_0,
        .filter_method = .method_0,
        .interlace_method = .none,
    };
    try std.testing.expectEqual(@as(IHDR.FromBytesError!IHDR.Data, data), IHDR.fromBytes(&data.toBytes()));
}

pub const PLTE = struct {
    pub const Palette = struct {
        entries_buf: [PLTE.max_entries]Entry,
        len_minus_one: LenMinusOne,

        pub fn AsBytes(comptime PalettePtr: type) type {
            var pointer_info = @typeInfo(PalettePtr).Pointer;
            assert(pointer_info.child == PLTE.Palette);
            assert(pointer_info.size == .One);

            pointer_info.child = u8;
            pointer_info.size = .Slice;
            return @Type(.{ .Pointer = pointer_info });
        }
        /// Returns a slice of bytes, backed by the palette entries.
        pub fn asBytes(palette: anytype) AsBytes(@TypeOf(palette)) {
            return std.mem.sliceAsBytes(entries(palette));
        }

        pub inline fn get(palette: *const PLTE.Palette, index: LenMinusOne) Entry {
            return palette.entries()[index];
        }

        pub fn Entries(comptime PalettePtr: type) type {
            var pointer_info = @typeInfo(PalettePtr).Pointer;
            assert(pointer_info.child == PLTE.Palette);
            assert(pointer_info.size == .One);

            pointer_info.child = Entry;
            pointer_info.size = .Slice;
            return @Type(.{ .Pointer = pointer_info });
        }
        pub inline fn entries(palette: anytype) Entries(@TypeOf(palette)) {
            return palette.entries_buf[0..len(palette.*)];
        }

        pub const Len = std.math.IntFittingRange(1, 256);
        pub inline fn len(data: PLTE.Palette) Len {
            return @as(Len, data.len_minus_one) + 1;
        }
    };

    pub const max_entries = 256;
    pub const Entry = extern struct { r: u8, g: u8, b: u8 };
    pub const LenMinusOne = std.math.IntFittingRange(1 - 1, 256 - 1);

    pub const FromBytesError = error{
        NoEntries,
        TooManyBytes,
        MissingOneByte,
        MissingTwoBytes,
    };
    pub fn fromBytes(bytes: []const u8) FromBytesError!PLTE.Palette {
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
        var result = PLTE.Palette{
            .entries_buf = undefined,
            .len_minus_one = @intCast(LenMinusOne, entries_slice.len - 1),
        };
        std.mem.copy(Entry, result.entries(), entries_slice);
        return result;
    }
};
test PLTE {
    // zig fmt: off
    try std.testing.expectError(error.NoEntries,       PLTE.fromBytes(&[_]u8{} ** 0));
    try std.testing.expectError(error.MissingTwoBytes, PLTE.fromBytes(&[_]u8{0} ** (256 * 3 - 2)));
    try std.testing.expectError(error.MissingOneByte,  PLTE.fromBytes(&[_]u8{0} ** (256 * 3 - 1)));
    try std.testing.expectError(error.TooManyBytes,    PLTE.fromBytes(&[_]u8{0} ** (256 * 3 + 1)));
    // zig fmt: on

    const data = PLTE.Palette{
        .entries_buf = [_]PLTE.Entry{.{ .r = 0, .g = 0, .b = 0 }} ** 256,
        .len_minus_one = 256 - 1,
    };
    try std.testing.expectEqual(@as(PLTE.FromBytesError!PLTE.Palette, data), PLTE.fromBytes(data.asBytes()));
}

pub const tRNS = struct {
    /// Data for `tRNS` when image is using color type 0
    pub const Grayscale = enum(u16) {
        _,

        pub const color_type = IHDR.ColorType.grayscale;

        pub inline fn fromBytes(bytes: [2]u8) Grayscale {
            const int = std.mem.readIntBig(WithBitDepth(.@"16"), &bytes);
            return fromInt(.@"16", int);
        }

        pub inline fn toBytes(gray_level: Grayscale) [2]u8 {
            var bytes: [2]u8 = undefined;
            std.mem.writeIntBig(u16, &bytes, @enumToInt(gray_level));
            return bytes;
        }

        pub inline fn fromInt(
            comptime bit_depth: IHDR.BitDepth,
            int: WithBitDepth(bit_depth),
        ) Grayscale {
            return @intToEnum(Grayscale, int);
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
        pub inline fn withBitDepth(
            transparency: Grayscale,
            comptime bit_depth: IHDR.BitDepth,
        ) ?WithBitDepth(bit_depth) {
            const int = @enumToInt(transparency);
            const result = std.math.cast(WithBitDepth(bit_depth), int) orelse return null;
            assert(int == result); // the discarded most significant bits should be 0
            return result;
        }
    };
    test Grayscale {
        const gray_level = tRNS.Grayscale.fromInt(.@"8", 255);
        try std.testing.expectEqual(gray_level, tRNS.Grayscale.fromBytes(gray_level.toBytes()));

        try std.testing.expectEqual(@as(?tRNS.Grayscale.WithBitDepth(.@"4"), null), gray_level.withBitDepth(.@"4"));
        // NOTE: Doesn't track the bit depth, so upcasting works incidentally.
        try std.testing.expectEqual(@as(?tRNS.Grayscale.WithBitDepth(.@"16"), 255), gray_level.withBitDepth(.@"16"));
    }

    /// Data for `tRNS` when image is using color type 2
    pub const Rgb = struct {
        r_level: u16,
        g_level: u16,
        b_level: u16,

        pub const color_type = IHDR.ColorType.rgb;

        pub inline fn fromBytes(bytes: [6]u8) Rgb {
            return Rgb{
                .r_level = std.mem.readIntBig(u16, bytes[0..2]),
                .g_level = std.mem.readIntBig(u16, bytes[2..4]),
                .b_level = std.mem.readIntBig(u16, bytes[4..6]),
            };
        }

        pub inline fn toBytes(transparency: Rgb) [6]u8 {
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
        pub inline fn withBitDepth(transparency: Rgb, comptime bit_depth: IHDR.BitDepth) ?WithBitDepth(bit_depth) {
            const FieldType = std.meta.FieldType;
            const Color = WithBitDepth(bit_depth);
            const result = Color{
                .r = std.math.cast(FieldType(Color, .r), transparency.r_level) orelse return null,
                .g = std.math.cast(FieldType(Color, .g), transparency.g_level) orelse return null,
                .b = std.math.cast(FieldType(Color, .b), transparency.b_level) orelse return null,
            };
            assert(result.r == transparency.r_level); // the discarded most significant bits should be 0
            assert(result.g == transparency.g_level); // the discarded most significant bits should be 0
            assert(result.b == transparency.b_level); // the discarded most significant bits should be 0
            return result;
        }
    };
    test Rgb {
        const rgb_level = tRNS.Rgb.fromBytes(std.mem.toBytes([3]u16{
            std.mem.nativeToBig(u16, 32),
            std.mem.nativeToBig(u16, 64),
            std.mem.nativeToBig(u16, 128),
        }));
        try std.testing.expectEqual(rgb_level, tRNS.Rgb.fromBytes(rgb_level.toBytes()));

        try std.testing.expectEqual(@as(?tRNS.Rgb.WithBitDepth(.@"16"), .{ .r = 32, .g = 64, .b = 128 }), rgb_level.withBitDepth(.@"16"));
        try std.testing.expectEqual(@as(?tRNS.Rgb.WithBitDepth(.@"2"), null), rgb_level.withBitDepth(.@"2"));
    }

    /// Data for `tRNS` when image is using color type 3
    pub const Palette = struct {
        entries_buf: [256]Entry,
        len_minus_one: LenMinusOne,

        pub const color_type = IHDR.ColorType.palette;

        pub const Entry = u8;
        pub const LenMinusOne = std.math.IntFittingRange(1 - 1, 256 - 1);

        pub const FromBytesError = error{
            NoEntries,
            TooManyBytes,
        };
        pub fn fromBytes(bytes: []const u8) FromBytesError!Palette {
            if (bytes.len == 0) return error.NoEntries;
            if (bytes.len > 256) return error.TooManyBytes;
            var result = Palette{
                .entries_buf = undefined,
                .len_minus_one = @intCast(LenMinusOne, bytes.len - 1),
            };
            std.mem.copy(Entry, result.entries(), bytes);
            return result;
        }

        pub fn AsBytes(comptime PalettePtr: type) type {
            var pointer_info = @typeInfo(PalettePtr).Pointer;
            assert(pointer_info.child == Palette);
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

        pub inline fn get(transparency: *const Palette, index: LenMinusOne) Entry {
            return transparency.entries()[index];
        }

        pub fn Entries(comptime PalettePtr: type) type {
            var pointer_info = @typeInfo(PalettePtr).Pointer;
            assert(pointer_info.child == Palette);
            assert(pointer_info.size == .One);

            pointer_info.child = Entry;
            pointer_info.size = .Slice;
            return @Type(.{ .Pointer = pointer_info });
        }
        pub inline fn entries(transparency: anytype) Entries(@TypeOf(transparency)) {
            return transparency.entries_buf[0..len(transparency.*)];
        }

        pub const Len = std.math.IntFittingRange(1, 256);
        pub inline fn len(data: Palette) Len {
            return @as(Len, data.len_minus_one) + 1;
        }
    };
    test Palette {
        try std.testing.expectError(error.NoEntries, tRNS.Palette.fromBytes(&[_]u8{}));
        try std.testing.expectError(error.TooManyBytes, tRNS.Palette.fromBytes(&[_]u8{0} ** 257));
        const palette = tRNS.Palette.fromBytes(&[_]u8{ 0, 1, 2, 3, 4, 5 }) catch undefined;
        _ = palette;
    }
};
comptime {
    if (builtin.is_test) {
        _ = tRNS;
    }
}

pub const TransparentPaletteIterator = struct {
    palette: *const PLTE.Palette,
    transparency: *const tRNS.Palette,
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
    palette: *const PLTE.Palette,
    transparency: *const tRNS.Palette,
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
    const transparency = try tRNS.Palette.fromBytes(&[_]u8{
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
pub const iCCP = struct {
    pub const profile_name_min_length = ReadNameResult.min_length;
    pub const profile_name_max_length = ReadNameResult.max_length;

    pub const ReadNameResult = shared.ReadBoundedStringResult(.{
        .min_length = 1,
        .max_length = 79,
    });

    /// Attempts to read the ICC profile name from the stream.
    ///
    /// The result of the function call is written to `out`,
    /// whether returning normally or with a write error.
    ///
    /// If the initialised tag is `ok`, the payload is the ICC profile name.
    /// Otherwise, the tag and payload describe the failure state.
    ///
    /// On success, consider calling `validateName*`, to confirm
    /// whether the result is a valid Latin-1 encoded ICC profile name.
    pub fn readName(
        out: *ReadNameResult,
        reader: anytype,
    ) @TypeOf(reader).Error!void {
        return @call(
            .always_inline,
            ReadNameResult.readString,
            .{ out, reader },
        );
    }

    pub const ValidateNameError = error{
        /// Encountered a byte representing a non-printable
        /// character, which also isn't a space.
        InvalidChar,
        /// Space encountered before any printable characters.
        LeadingSpace,
        /// Space encountered immediately following another space.
        ConsecutiveSpace,
        /// Space encountered at the end of the string.
        TrailingSpace,
    };
    /// Validates a slice of bytes returned by `readName`.
    pub fn validateName(
        /// Should be from `ReadNameResult.ok`, which was written by a call to `readName`.
        name: []const u8,
    ) ValidateNameError!void {
        var index = {};
        return @call(.always_inline, validateNameTemplate, .{ name, false, &index });
    }
    /// Like `validateName`, but takes an out parameter `index`.
    pub fn validateNameIndex(
        /// Should be from `ReadNameResult.ok`, which was written by a call to `readName`.
        name: []const u8,
        /// The index in `name` at which the error occurred, if an error occured.
        index: *usize,
    ) ValidateNameError!void {
        return @call(.always_inline, validateNameTemplate, .{ name, true, index });
    }
    fn validateNameTemplate(
        name: []const u8,
        comptime want_index: bool,
        index: *if (want_index) usize else void,
    ) ValidateNameError!void {
        assert(name.len != 0); // shouldn't be possible if this is from `readName`
        assert(name.len <= profile_name_max_length); // shoudln't be possible if this is from `readName`
        index.* = undefined;

        @as(error{ LeadingSpace, InvalidChar }!void, switch (name[0]) {
            32 => error.LeadingSpace,
            33...126,
            161...255,
            => {},
            else => error.InvalidChar,
        }) catch |err| {
            index.* = if (want_index) 0 else {};
            return err;
        };

        if (name.len == 1) return;

        @as(error{ TrailingSpace, InvalidChar }!void, switch (name[name.len - 1]) {
            32 => error.TrailingSpace,
            33...126,
            161...255,
            => {},
            else => error.InvalidChar,
        }) catch |err| {
            index.* = if (want_index) name.len - 1 else {};
            return err;
        };

        var state: enum { print, space } = .print;

        for (name[1..]) |char, i| {
            errdefer index.* = if (want_index) i + 1 else {};
            switch (state) {
                .print => switch (char) {
                    32 => state = .space,
                    33...126,
                    161...255,
                    => {},
                    else => return error.InvalidChar,
                },
                .space => switch (char) {
                    32 => return error.ConsecutiveSpace,
                    33...126,
                    161...255,
                    => state = .print,
                    else => return error.InvalidChar,
                },
            }
        }
    }

    pub const CompressionMethod = enum(u8) {
        /// zlib datastream with deflate compression
        method_0,

        pub inline fn toByte(method: CompressionMethod) u8 {
            return @enumToInt(method);
        }

        pub inline fn fromByte(byte: u8) ?CompressionMethod {
            return util.intToEnum(CompressionMethod, byte);
        }
    };
};
test iCCP {
    var fbs = std.io.fixedBufferStream(&[_]u8{ 'f', 'o', 'o', ' ', 'b', 'a', 'r', 0, @enumToInt(iCCP.CompressionMethod.method_0) });
    const fbr = fbs.reader();

    var maybe_name: iCCP.ReadNameResult = undefined;
    iCCP.readName(&maybe_name, fbr) catch |err| switch (err) {};
    try std.testing.expectEqual(iCCP.ReadNameResult.ok, maybe_name);
    try std.testing.expectEqual(@as(iCCP.ValidateNameError!void, {}), iCCP.validateName(maybe_name.ok.constSlice()));
    try std.testing.expectEqualStrings(maybe_name.ok.constSlice(), "foo bar");

    try std.testing.expectEqual(@as(?iCCP.CompressionMethod, .method_0), iCCP.CompressionMethod.fromByte(try fbr.readByte()));

    var index: usize = undefined;

    try std.testing.expectError(error.InvalidChar, iCCP.validateNameIndex(&[_]u8{ 'f', 'o', 'o', ' ', 'b', 0, 'r' }, &index));
    try std.testing.expectEqual(@as(usize, 5), index);

    try std.testing.expectError(error.ConsecutiveSpace, iCCP.validateNameIndex(&[_]u8{ 'f', 'o', 'o', ' ', ' ', 'b', 'r' }, &index));
    try std.testing.expectEqual(@as(usize, 4), index);

    try std.testing.expectError(error.LeadingSpace, iCCP.validateNameIndex(&[_]u8{ ' ', ' ', ' ', 'b' }, &index));
    try std.testing.expectEqual(@as(usize, 0), index);

    try std.testing.expectError(error.TrailingSpace, iCCP.validateNameIndex(&[_]u8{ 'b', 'o', 'c', 'a', 'd', 'i', ' ' }, &index));
    try std.testing.expectEqual(@as(usize, 6), index);
}

/// Returns a formatter for a Latin-1 string.
/// This may be useful, because Latin-1 is encoding
/// of the textual data in various PNG chunks.
pub inline fn fmtLatin1(latin1_string: []const u8) FmtLatin1 {
    return .{ .latin1_string = latin1_string };
}
pub const FmtLatin1 = struct {
    latin1_string: []const u8,

    pub fn format(
        formatter: FmtLatin1,
        comptime fmt_str: []const u8,
        options: std.fmt.FormatOptions,
        writer: anytype,
    ) @TypeOf(writer).Error!void {
        comptime assert(fmt_str.len == 0);
        _ = options;
        for (formatter.latin1_string) |char| {
            var cp_buf: [4]u8 = undefined;
            // none of the errors returned by this function should
            // be applicable to codepoints this small,
            // so mark as unreachable.
            const cp_byte_len = std.unicode.utf8Encode(char, &cp_buf) catch unreachable;
            try writer.writeAll(cp_buf[0..cp_byte_len]);
        }
    }
};
test fmtLatin1 {
    try std.testing.expectFmt("abc", "{}", .{fmtLatin1("abc")}); // ascii strings obviously work
    try std.testing.expectFmt("ß Ø £", "{}", .{fmtLatin1(&[_]u8{ 'ß', ' ', 'Ø', ' ', '£' })});
}

const shared = struct {
    const ReadBoundedStringResultTag = enum {
        ok,
        bad_sentinel,
        no_sentinel,
        incomplete,
        empty_string,
        empty_stream,
        immediate_error,
    };
    const ReadBoundedStringResultConfig = struct {
        min_length: comptime_int,
        max_length: comptime_int,
    };
    fn ReadBoundedStringResult(
        comptime config: ReadBoundedStringResultConfig,
    ) type {
        return union(ReadBoundedStringResultTag) {
            const Result = @This();
            /// the stream returned the whole string, including the null
            /// byte separator (not included in the attached result).
            ok: NameBytes,
            /// the stream returned the maximum number of bytes allowed
            /// for the profile name, but returned a non-zero byte when
            /// the null byte separator was expected.
            bad_sentinel: BadSentinel,
            /// the stream returned the maximum number of bytes allowed
            /// for the profile name, but didn't return a null byte (either
            /// the stream returned an error or it ended).
            /// Attached is the array of bytes that the stream returned beforehand.
            no_sentinel: NoSentinel,
            /// the stream ended or returned an error before returning
            /// the null byte separator. Attached is a bounded array of
            /// what the stream returned beforehand.
            incomplete: Incomplete,
            /// the stream immedately returned the null byte separator,
            /// yielding an empty name.
            empty_string,
            /// the stream immedately ended.
            empty_stream,
            /// the stream immediately returned an error.
            immediate_error,

            const min_length = config.min_length;
            const max_length = config.max_length;

            pub const NameBytes = std.BoundedArray(u8, max_length);
            pub const BadSentinel = struct {
                name: [max_length]u8,
                sentinel: u8,
            };
            pub const NoSentinel = [max_length]u8;
            pub const Incomplete = std.BoundedArray(u8, max_length - 1);

            // -- convenience --

            pub const UnwrapError = error{
                BadSentinel,
                NoSentinel,
                Incomplete,
                EmptyString,
                EmptyStream,
                ImmediateError,
            };
            /// Converts this union to an error union.
            /// Tags other than `ok` are errors.
            pub fn unwrap(result: Result) UnwrapError!NameBytes {
                return switch (result) {
                    .ok => |name| name,
                    .bad_sentinel => error.BadSentinel,
                    .no_sentinel => error.NoSentinel,
                    .incomplete => error.Incomplete,
                    .empty_string => error.EmptyString,
                    .empty_stream => error.EmptyStream,
                    .immediate_error => error.ImmediateError,
                };
            }

            /// Attempts to read a bounded string from the stream,
            /// which is expected to be followed by a null byte.
            ///
            /// The result of the function call is written to `out`,
            /// whether returning normally or with a write error.
            ///
            /// If the initialised tag is `ok`, the payload is the ICC profile name.
            /// Otherwise, the tag and payload describe the failure state.
            fn readString(
                out: *Result,
                reader: anytype,
            ) @TypeOf(reader).Error!void {
                out.* = undefined;
                // just to check the result is a valid value
                // by the end of the function call
                defer switch (out.*) {
                    else => {},
                };

                var name: Result.NameBytes = .{};
                if (util.readByteOrNull(reader) catch |err| {
                    out.* = .immediate_error;
                    return err;
                }) |first_byte| {
                    if (first_byte == 0) {
                        out.* = .empty_string;
                        return;
                    }
                    name.appendAssumeCapacity(first_byte);
                } else {
                    out.* = .empty_stream;
                    return;
                }

                while (name.len < max_length) {
                    const maybe_byte = util.readByteOrNull(reader) catch |err| {
                        out.* = .{ .incomplete = Result.Incomplete.fromSlice(name.constSlice()) catch unreachable };
                        return err;
                    };
                    const byte: u8 = maybe_byte orelse {
                        out.* = .{ .incomplete = Result.Incomplete.fromSlice(name.constSlice()) catch unreachable };
                        return;
                    };
                    if (byte == 0) break;
                    name.appendAssumeCapacity(byte);
                } else {
                    assert(name.len == max_length);
                    const maybe_byte = util.readByteOrNull(reader) catch |err| {
                        out.* = .{ .no_sentinel = name };
                        return err;
                    };
                    const sentinel_byte: u8 = maybe_byte orelse {
                        out.* = .{ .no_sentinel = name.constSlice()[0..max_length].* };
                        return;
                    };
                    if (sentinel_byte != 0) {
                        out.* = .{ .bad_sentinel = Result.BadSentinel{
                            .name = name.constSlice()[0..max_length].*,
                            .sentinel = sentinel_byte,
                        } };
                        return;
                    }
                }

                out.* = .{ .ok = name };
                return;
            }
        };
    }
};

const util = struct {
    /// Read a byte from the stream, or return null if the stream is empty.
    inline fn readByteOrNull(reader: anytype) @TypeOf(reader).Error!?u8 {
        var byte: [1]u8 = undefined;
        return switch (try reader.read(&byte)) {
            0 => null,
            1 => @bitCast(u8, byte),
            else => unreachable,
        };
    }

    fn intToEnum(comptime Enum: type, int: @typeInfo(Enum).Enum.tag_type) ?Enum {
        inline for (@typeInfo(Enum).Enum.fields) |field| {
            if (int == field.value)
                return @field(Enum, field.name);
        }
        return null;
    }
};
