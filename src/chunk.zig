const std = @import("std");
const builtin = @import("builtin");
const assert = std.debug.assert;

const png = @import("main.zig");
const chunk = @This();

pub const Header = struct {
    length: u31,
    type: Type,

    pub const ValidateError = error{
        InvalidLength,
        InvalidType,
    };
    pub fn fromBytes(bytes: *const [8]u8) Raw {
        var result: Raw = @bitCast(Raw, bytes.*);
        result.length = std.mem.bigToNative(u32, result.length);
        return result;
    }

    pub fn toRaw(chunk_header: Header) Raw {
        return Raw{
            .length = chunk_header.length,
            .type = chunk_header.type.string(),
        };
    }

    pub const Raw = extern struct {
        length: u32,
        type: [4]u8,

        pub fn validate(raw: Raw) ValidateError!Header {
            return Header{
                .length = try (std.math.cast(u31, raw.length) orelse error.InvalidLength),
                .type = try (Type.from(raw.type) orelse error.InvalidType),
            };
        }

        pub fn toBytes(raw: Raw) [8]u8 {
            var copy = raw;
            copy.length = std.mem.nativeToBig(u32, copy.length);
            return @bitCast([8]u8, copy);
        }
    };
};

test Header {
    const header = Header{
        .length = 13,
        .type = .IHDR,
    };
    try std.testing.expectEqual(@as(Header.ValidateError!Header, header), Header.fromBytes(&header.toRaw().toBytes()).validate());
}

pub const Type = enum(u32) {
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

    pub inline fn from(bytes: [4]u8) ?Type {
        const ct = @intToEnum(Type, valueFromString(bytes));
        if (!ct.isValidAscii()) return null;
        return ct;
    }

    pub inline fn string(ct: Type) [4]u8 {
        return @bitCast([4]u8, @enumToInt(ct));
    }

    pub inline fn isValidAscii(ct: Type) bool {
        const bytes: [4]u8 = @bitCast([4]u8, @enumToInt(ct));
        inline for (bytes) |byte| {
            if (!std.ascii.isAlphabetic(byte)) return false;
        }
        return true;
    }

    pub inline fn isAncillary(ct: Type) bool {
        return propertyBitMask(ct, 0) != 0;
    }
    pub inline fn isPrivate(ct: Type) bool {
        return propertyBitMask(ct, 1) != 0;
    }
    pub fn isReserved(ct: Type) bool {
        return propertyBitMask(ct, 2) != 0;
    }
    pub inline fn isSafeToCopy(ct: Type) bool {
        return propertyBitMask(ct, 3) != 0;
    }

    inline fn propertyBitMask(ct: Type, index: u2) u8 {
        return string(ct)[index] & 0b00100000;
    }

    inline fn valueFromString(str: [4]u8) u32 {
        return @bitCast(u32, str);
    }
};

test Type {
    try std.testing.expectEqual(@as(?Type, null), Type.from("0ABC".*));
    try std.testing.expect(Type.from("fOOb".*) != null);
    try std.testing.expectEqual(@as(?Type, Type.IHDR), Type.from("IHDR".*));
    inline for (@typeInfo(Type).Enum.fields) |field| {
        try std.testing.expect(field.name.len == 4);
        try std.testing.expect(field.value == @bitCast(u32, field.name[0..4].*));

        const tag: Type = @intToEnum(Type, field.value);
        try std.testing.expect(Type.isValidAscii(tag));
        try std.testing.expect(!Type.isReserved(tag)); // only support non-reserved chunks among the named fields, same as the target PNG version
        try std.testing.expect(!Type.isPrivate(tag)); // only allow for public chunk types among the named fields
    }
}

/// Returns an iterator that reads bytes from the given stream,
/// and outputs chunk headers, providing a wrapped stream that
/// the caller should use to retrieve the bytes following each
/// chunk (the stream will make sure it only returns the content
/// of the current chunk).
pub inline fn iterator(
    reader: anytype,
    /// Enables calculating the CRC code.
    comptime calculate_crc: bool,
) Iterator(@TypeOf(reader), calculate_crc) {
    return .{ .inner = reader };
}
pub const IteratorSignatureCheck = union(enum) {
    ok,
    bad_signature: [png.signature.len]u8,
    incomplete_signature: std.BoundedArray(u8, png.signature.len - 1),

    pub const UnwrapError = error{ BadPngSignature, IncompletePngSignature };
    pub inline fn unwrap(result: IteratorSignatureCheck) UnwrapError!void {
        return switch (result) {
            .ok => {},
            .bad_signature => error.BadPngSignature,
            .incomplete_signature => error.IncompletePngSignature,
        };
    }
};
pub const IteratorNextResult = union(enum) {
    none,
    ok: Header,
    incomplete_chunk_header: std.BoundedArray(u8, 8 - 1),
    invalid_chunk_header: Header.Raw,

    pub const UnwrapError = error{ IncompleteBytes, InvalidChunkHeader };
    pub inline fn unwrap(result: IteratorNextResult) UnwrapError!?Header {
        return switch (result) {
            .none => null,
            .ok => |header| header,
            .incomplete_chunk_header => error.IncompleteBytes,
            .invalid_chunk_header => error.InvalidChunkHeader,
        };
    }
};
pub fn Iterator(
    comptime ReaderType: type,
    /// Enables calculating the CRC code.
    comptime calculate_crc: bool,
) type {
    return struct {
        const Self = @This();
        inner: Inner,
        remaining_bytes: u31 = undefined,
        crc_hasher: if (calculate_crc) std.hash.Crc32 else void = if (calculate_crc) std.hash.Crc32.init() else {},
        expected_crc: u32 = undefined,
        state: State = .start,

        pub const Inner = ReaderType;

        pub const SignatureCheck = IteratorSignatureCheck;
        /// The result describes whether the PNG
        /// signature was attained correctly from
        /// the stream.
        ///
        /// On successful return, the caller can call `next*`.
        pub fn checkSignature(iter: *Self) Inner.Error!SignatureCheck {
            var out: SignatureCheck = undefined;
            try @call(.always_inline, checkSignatureAdvanced, .{ iter, &out });
            return out;
        }
        /// The result is written to `out`, whether
        /// returning normally or with an error. It
        /// describes whether the PNG signature was
        /// attained correctly from the stream.
        ///
        /// On successful return, the caller can call `next*`.
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

            var bytes: std.BoundedArray(u8, png.signature.len) = .{};
            iter.inner.readIntoBoundedBytes(png.signature.len, &bytes) catch |err| {
                assert(bytes.len < png.signature.len);
                out.* = .{ .incomplete_signature = std.BoundedArray(u8, png.signature.len).fromSlice(bytes.constSlice()) catch unreachable };
                iter.state = .end;
                return err;
            };
            if (bytes.len < png.signature.len) {
                out.* = .{ .incomplete_signature = std.BoundedArray(u8, png.signature.len - 1).fromSlice(bytes.constSlice()) catch unreachable };
                iter.state = .end;
                return;
            }
            if (!std.mem.eql(u8, bytes.constSlice(), &png.signature)) {
                out.* = .{ .bad_signature = bytes.constSlice()[0..png.signature.len].* };
                iter.state = .end;
                return;
            }

            out.* = .ok;
            iter.state = .expecting_header;
            return;
        }

        pub const NextResult = IteratorNextResult;
        /// The result describes whether the chunk
        /// header was attained from the stream and
        /// decoded correctly.
        ///
        /// On successful return, the caller should exhaust
        /// the stream returned by `dataReader`.
        pub fn next(iter: *Self) Inner.Error!NextResult {
            var out: NextResult = undefined;
            try @call(.always_inline, nextAdvanced, .{ iter, &out });
            return out;
        }
        /// The result is written to `out`, whether
        /// returning normally or with an error. It
        /// describes whether the chunk header was
        /// attained from the stream and decoded correctly.
        ///
        /// On successful return, the caller should exhaust
        /// the stream returned by `dataReader`.
        pub fn nextAdvanced(iter: *Self, out: *NextResult) Inner.Error!void {
            switch (iter.state) {
                .start => unreachable, // meant to call `checkSignature*` first.s
                .end => {
                    out.* = .none;
                    return;
                },
                .expecting_header => {},
                .exhausting_data => unreachable, // supposed to call meant to call `fetchExpectedCrc` first.
                .got_expected_crc => {
                    assert(iter.remaining_bytes == 0);
                    iter.remaining_bytes = undefined;
                    iter.expected_crc = undefined;
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
            const raw_reader = Header.fromBytes(bytes.constSlice()[0..8]);
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

        pub const DataReader = std.io.Reader(*Self, Inner.Error, readData);
        /// Returns a stream that returns the raw contents
        /// of the current chunk.
        ///
        /// The caller should exhaust this stream, and
        /// then call `fetchExpectedCrc`.
        /// The caller may also call `getActualCrc`
        /// after exhausting this stream.
        pub inline fn dataReader(iter: *Self) DataReader {
            return .{ .context = iter };
        }

        /// Returns true on successfully reading the CRC code,
        /// false if the stream ended before returning it fully.
        /// Returning a read error implies the same as returning false.
        ///
        /// The caller may also call `getExpectedCrc`
        /// after calling this function.
        /// The caller may also call `next*` after
        /// calling this function, but not before `getExpectedCrc`.
        pub fn fetchExpectedCrc(iter: *Self) Inner.Error!bool {
            switch (iter.state) {
                .start,
                .end,
                .expecting_header,
                => unreachable, //  not valid to call yet

                .got_expected_crc => unreachable, // shouldn't call twice
                .exhausting_data => {
                    assert(iter.remaining_bytes == 0); // caller is meant to exhaust the stream returned by `.dataReader()`.
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

        /// Get the CRC code which was calculated
        /// based off of the data read from the stream.
        pub inline fn getActualCrc(iter: *const Self) u32 {
            comptime assert(calculate_crc); // must enable CRC calculation in order to use this
            switch (iter.state) {
                .start,
                .end,
                .expecting_header,
                => unreachable, //  not valid to call yet

                .exhausting_data,
                .got_expected_crc,
                => {},
            }
            assert(iter.remaining_bytes == 0); // must exhaust the data reader stream first
            var hasher_copy: std.hash.Crc32 = iter.crc_hasher;
            return hasher_copy.final();
        }

        /// Get the CRC code found after the chunk data.
        /// Must only be called after calling `fetchExpectedCrc`.
        pub inline fn getExpectedCrc(iter: *const Self) u32 {
            switch (iter.state) {
                .start,
                .end,
                .expecting_header,
                => unreachable, //  not valid to call yet

                .exhausting_data => unreachable, // must first call `fetchExpectedCrc`.
                .got_expected_crc => {},
            }
            assert(iter.remaining_bytes == 0); // the caller probably forgot to exhaust the `dataReader` stream.
            return iter.expected_crc;
        }

        const State = enum {
            start,
            end,
            expecting_header,
            exhausting_data,
            got_expected_crc,
        };

        fn readData(iter: *Self, buffer: []u8) Inner.Error!usize {
            switch (iter.state) {
                .start,
                .end,
                .expecting_header,
                .got_expected_crc,
                => unreachable, // should only read from this stream after calling `next`.
                .exhausting_data => {},
            }
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
    pub const Entry = extern struct { r: u8, g: u8, b: u8 };
    comptime {
        assert(@sizeOf(Entry) == 3);
        assert(@offsetOf(Entry, "r") == 0);
        assert(@offsetOf(Entry, "g") == 1);
        assert(@offsetOf(Entry, "b") == 2);
    }

    pub const min_entries = 1;
    pub const max_entries = 256;

    pub const BytesAsEntriesError = error{
        NoEntries,
        TooManyBytes,
        MissingOneByte,
        MissingTwoBytes,
    };
    pub fn bytesAsEntries(bytes: []const u8) BytesAsEntriesError![]const Entry {
        if (bytes.len == 0) return error.NoEntries;
        if (bytes.len > 256 * 3) return error.TooManyBytes;
        try switch (bytes.len % 3) {
            0 => {},
            1 => error.MissingTwoBytes,
            2 => error.MissingOneByte,
            else => unreachable,
        };
        return std.mem.bytesAsSlice(Entry, bytes);
    }
};
test PLTE {
    // zig fmt: off
    try std.testing.expectError(error.NoEntries,       PLTE.bytesAsEntries(&[_]u8{ } ** 0));
    try std.testing.expectError(error.MissingTwoBytes, PLTE.bytesAsEntries(&[_]u8{0} ** (256 * 3 - 2)));
    try std.testing.expectError(error.MissingOneByte,  PLTE.bytesAsEntries(&[_]u8{0} ** (256 * 3 - 1)));
    try std.testing.expectError(error.TooManyBytes,    PLTE.bytesAsEntries(&[_]u8{0} ** (256 * 3 + 1)));
    // zig fmt: on
    try std.testing.expectEqualSlices(
        PLTE.Entry,
        &[_]PLTE.Entry{ .{ .r = 1, .g = 2, .b = 3 }, .{ .r = 22, .g = 0, .b = 44 } },
        PLTE.bytesAsEntries(&[_]u8{ 1, 2, 3 } ++ [_]u8{ 22, 0, 44 }) catch @panic("unreachable"),
    );
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
        /// Red:   2 bytes, range 0 .. (2^bitdepth)-1
        r_level: u16,
        /// Green: 2 bytes, range 0 .. (2^bitdepth)-1
        g_level: u16,
        /// Blue:  2 bytes, range 0 .. (2^bitdepth)-1
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

    pub const color_type_3 = struct {
        pub const min_entries = 1;
        pub const max_entries = 256;

        pub const RgbaEntry = struct {
            r: u8,
            g: u8,
            b: u8,
            a: u8,
        };
        pub const PaletteIterator = struct {
            rgb_palette: []const PLTE.Entry,
            alpha_palette: []const u8,
            index: u8 = 0,

            pub fn next(iter: *PaletteIterator) ?RgbaEntry {
                if (iter.index == iter.rgb_palette.len) return null;
                const rgb = iter.rgb_palette[iter.index];
                var result = RgbaEntry{
                    .r = rgb.r,
                    .g = rgb.g,
                    .b = rgb.b,
                    .a = 255,
                };
                if (iter.index < iter.alpha_palette.len) {
                    result.a = iter.alpha_palette[iter.index];
                }
                iter.index += 1;
                return result;
            }
        };

        /// Iterates over the entries of a `PLTE` and a `tRNS`,
        /// outputting combined RGBA entries.
        pub fn paletteIterator(
            rgb_palette: []const PLTE.Entry,
            alpha_palette: []const u8,
        ) error{TooManyTransparencyEntries}!PaletteIterator {
            assert(rgb_palette.len >= PLTE.min_entries);
            assert(rgb_palette.len <= PLTE.max_entries);

            assert(alpha_palette.len >= color_type_3.min_entries);
            assert(alpha_palette.len <= color_type_3.max_entries);

            if (rgb_palette.len < alpha_palette.len)
                return error.TooManyTransparencyEntries;

            return PaletteIterator{
                .rgb_palette = rgb_palette,
                .alpha_palette = alpha_palette,
            };
        }
    };

    test color_type_3 {
        const palette: []const PLTE.Entry = try PLTE.bytesAsEntries(&[_]u8{
            0x00, 0x00, 0xff,
            0xff, 0xff, 0x00,
            0x32, 0x32, 0x32,
            0xff, 0x00, 0x0,
        });
        const transparency: []const u8 = &[_]u8{
            0xef,
            0xab,
            0xfa,
        };
        var iter = try color_type_3.paletteIterator(palette, transparency);
        const Expected = ?color_type_3.RgbaEntry;
        try std.testing.expectEqual(@as(Expected, .{ .r = 0x00, .g = 0x00, .b = 0xff, .a = 0xef }), iter.next());
        try std.testing.expectEqual(@as(Expected, .{ .r = 0xff, .g = 0xff, .b = 0x00, .a = 0xab }), iter.next());
        try std.testing.expectEqual(@as(Expected, .{ .r = 0x32, .g = 0x32, .b = 0x32, .a = 0xfa }), iter.next());
        try std.testing.expectEqual(@as(Expected, .{ .r = 0xff, .g = 0x00, .b = 0x00, .a = 0xff }), iter.next());
        try std.testing.expectEqual(@as(Expected, null), iter.next());
        try std.testing.expectEqual(@as(Expected, null), iter.next());
    }
};
comptime {
    if (builtin.is_test) {
        _ = tRNS;
    }
}

pub const gAMA = struct {
    //! Gamma: 4 bytes

    /// A typed integer representing the value of the `gAMA` chunk.
    pub const Value = enum(u32) {
        _,

        pub inline fn toBytes(gamma: gAMA.Value) [4]u8 {
            var bytes: [4]u8 = undefined;
            const int: u32 = @enumToInt(gamma);
            std.mem.writeIntBig(u32, &bytes, int);
            return bytes;
        }

        pub inline fn times_100_000(gamma: gAMA.Value) u32 {
            return @enumToInt(gamma);
        }

        /// The actual gamma value in floating point representation.
        pub inline fn as(gamma: gAMA.Value, comptime T: type) T {
            return @intToFloat(T, gamma.times_100_000()) / 100_000.0;
        }
    };

    pub inline fn fromBytes(bytes: [4]u8) gAMA.Value {
        const int = std.mem.readIntBig(u32, &bytes);
        return gAMA.fromInt(int);
    }
    pub inline fn fromInt(int_times_100k: u32) gAMA.Value {
        return @intToEnum(gAMA.Value, int_times_100k);
    }
};
test gAMA {
    const gamma: gAMA.Value = sRGB.compat_gAMA();
    try std.testing.expectEqual(gamma, gAMA.fromBytes(gamma.toBytes()));
    try std.testing.expectEqual(@as(f128, 0.45_455), gamma.as(f128));
    try std.testing.expectEqual(@as(u32, 45_455), gamma.times_100_000());
}

pub const cHRM = struct {
    //! White Point x: 4 bytes
    //! White Point y: 4 bytes
    //! Red x:         4 bytes
    //! Red y:         4 bytes
    //! Green x:       4 bytes
    //! Green y:       4 bytes
    //! Blue x:        4 bytes
    //! Blue y:        4 bytes

    pub const Data = struct {
        white_point_x: ValueTimes100K,
        white_point_y: ValueTimes100K,
        red_x: ValueTimes100K,
        red_y: ValueTimes100K,
        green_x: ValueTimes100K,
        green_y: ValueTimes100K,
        blue_x: ValueTimes100K,
        blue_y: ValueTimes100K,

        pub fn toBytes(chroma: cHRM.Data) [32]u8 {
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
    };

    pub fn fromBytes(bytes: *const [32]u8) cHRM.Data {
        return cHRM.Data{
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
    const chroma: cHRM.Data = sRGB.compat_cHRM();
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

pub const sRGB = struct {
    /// Represents the rendering intent of the sRGB chunk.
    pub const RenderingIntent = enum(u8) {
        perceptual = 0,
        relative_colorimetric = 1,
        saturation = 2,
        absolute_colorimetric = 3,

        pub inline fn toByte(render_intent: RenderingIntent) u8 {
            return @enumToInt(render_intent);
        }
    };

    pub inline fn fromByte(byte: u8) ?sRGB.RenderingIntent {
        return util.intToEnum(sRGB.RenderingIntent, byte);
    }

    /// The only `gAMA` chunk value allowed for an application writing both
    /// an `sRGB` chunk and a `gAMA` chunk for compatibility.
    pub inline fn compat_gAMA() gAMA.Value {
        comptime return gAMA.fromInt(45455);
    }

    /// The only `cHRM` chunk value allowed for an application writing both
    /// an `sRGB` chunk and a `cHRM` chunk for compatibility.
    pub inline fn compat_cHRM() cHRM.Data {
        comptime return cHRM.Data{
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
    const color_space: sRGB.RenderingIntent = .absolute_colorimetric;
    try std.testing.expectEqual(@as(?sRGB.RenderingIntent, color_space), sRGB.fromByte(color_space.toByte()));
}

pub const iCCP = struct {
    //! Profile name:       1-79 bytes (character string)
    //! Null separator:     1 byte
    //! Compression method: 1 byte
    //! Compressed profile: n bytes

    pub const name_min_length = 1;
    pub const name_max_length = 79;

    pub const ReadNameResult = shared.ReadBoundedStringResult(.{
        .min_length = name_min_length,
        .max_length = name_max_length,
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
        assert(name.len <= name_max_length); // shoudln't be possible if this is from `readName`
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

pub const tEXt = struct {
    //! Keyword:        1-79 bytes (character string)
    //! Null separator: 1 byte
    //! Text:           n bytes (character string)

    pub const keyword_min_length = 1;
    pub const keyword_max_length = 79;

    pub const ReadKeywordResult = shared.ReadBoundedStringResult(.{
        .min_length = keyword_min_length,
        .max_length = keyword_max_length,
    });

    /// Attempts to read the keyword of the tEXt chunk.
    ///
    /// The result of the function call is written to `out`,
    /// whether returning normally or with a write error.
    ///
    /// If the initialised tag is `ok`, the payload is the keyword.
    /// Otherwise, the tag and payload describe the failure state.
    pub fn readKeyword(
        out: *ReadKeywordResult,
        reader: anytype,
    ) @TypeOf(reader).Error!void {
        return @call(
            .always_inline,
            ReadKeywordResult.readString,
            .{ out, reader },
        );
    }
};
test tEXt {
    _ = tEXt.ReadKeywordResult;
    _ = tEXt.readKeyword;
    return error.SkipZigTest;
}

/// Returns a formatter for a Latin-1 string,
/// to print it in UTF-8 encoding.
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
