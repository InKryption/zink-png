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
        _ = fmt_str;
        _ = fmt_options;
        const lazy = struct {
            const unknown_format_str = @compileError("Unknown format string '" ++ fmt_str ++ "'.\n");
        };

        if (fmt_str.len != 0) lazy.unknown_format_str;
        return writer.print(@typeName(@This()) ++ "({s})", .{ch_type.string()});
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

    pub fn fromBytes(bytes: *const [8]u8) ChunkHeader {
        return switch (ChunkHeader.parseBuffer(bytes)) {
            .ok => |value| value,
            else => unreachable,
        };
    }

    pub const ParseBufferResult = union(enum) {
        ok: ChunkHeader,
        no_type: NoType,
        invalid_length: InvalidLength,
        no_length: NoLength,

        pub const NoType = struct { length: u31 };
        pub const InvalidLength = struct { length: u32 };
        pub const NoLength = struct {};
    };

    pub fn parseBuffer(buffer: []const u8) ParseBufferResult {
        var fbs = std.io.fixedBufferStream(buffer);
        return switch (ChunkHeader.parseReader(fbs.reader())) {
            .ok => |value| ParseBufferResult{ .ok = value },
            .no_type_eos => |info| ParseBufferResult{ .no_type = info },
            .no_type_err => unreachable,
            .invalid_length => |info| ParseBufferResult{ .invalid_length = info },
            .no_length_eos => |info| ParseBufferResult{ .no_length = info },
            .no_length_err => unreachable,
        };
    }

    pub const ParseReaderResultTag = enum {
        ok,
        no_type_eos,
        no_type_err,
        invalid_length,
        no_length_eos,
        no_length_err,
    };
    pub fn ParseReaderResult(comptime ReaderError: type) type {
        return union(ParseReaderResultTag) {
            /// success
            ok: ChunkHeader,
            /// encountered end of stream while trying to read type
            no_type_eos: NoTypeEos,
            /// encountered error while trying to read type
            no_type_err: NoTypeErr,
            /// encountered an invalid length
            invalid_length: InvalidLength,
            /// encountered end of stream while trying to read length
            no_length_eos: NoLengthEos,
            /// encountered error while trying to read length
            no_length_err: NoLengthErr,

            pub const ReadError = ReaderError;
            pub const NoTypeEos = ParseBufferResult.NoType;
            pub const NoTypeErr = struct { length: u31, err: ReadError };
            pub const InvalidLength = ParseBufferResult.InvalidLength;
            pub const NoLengthEos = ParseBufferResult.NoLength;
            pub const NoLengthErr = struct { err: ReadError };
        };
    }

    /// Returns null if the stream ends before returning the required number of bytes for a chunk header.
    pub fn parseReader(reader: anytype) ParseReaderResult(util.NormalizedErrorSet(@TypeOf(reader).Error)) {
        const PResult = ParseReaderResult(util.NormalizedErrorSet(@TypeOf(reader).Error));

        const length = if (util.io.readIntBigOrNull(reader, u32)) |maybe_length| length: {
            const length = maybe_length orelse return PResult{ .no_length_eos = .{} };
            break :length std.math.cast(u31, length) orelse
                return PResult{ .invalid_length = .{ .length = length } };
        } else |err| {
            return PResult{ .no_length_err = .{ .err = err } };
        };

        const type_value = util.io.readIntBigOrNull(reader, u32) catch |err| {
            return PResult{ .no_type_err = .{ .err = err, .length = length } };
        } orelse return PResult{ .no_type_eos = .{ .length = length } };

        return PResult{ .ok = ChunkHeader{
            .length = length,
            .type = @intToEnum(ChunkType, type_value),
        } };
    }
};

test "ChunkType & ChunkHeader" {
    try std.testing.expectEqualStrings("WOAH", &ChunkType.from("WOAH").string());
    for (@as([4096]void, undefined)) |_, i| {
        inline for (comptime std.enums.values(ChunkType)) |ch_type| {
            const header = ChunkHeader{
                .length = @intCast(u31, i),
                .type = ch_type,
            };
            try std.testing.expectEqual(header, ChunkHeader.fromBytes(&header.toBytes()));
        }
    }

    var fbs: std.io.FixedBufferStream([]const u8) = undefined;
    var elr: util.io.ErrorLimitedReader(@TypeOf(fbs).Reader) = undefined;

    fbs = std.io.fixedBufferStream("");
    elr = undefined;
    try std.testing.expectEqual(ChunkHeader.ParseReaderResultTag.no_length_eos, ChunkHeader.parseReader(fbs.reader()));

    fbs = std.io.fixedBufferStream("");
    elr = util.io.errorLimitedReader(fbs.reader(), 0);
    try std.testing.expectEqual(ChunkHeader.ParseReaderResultTag.no_length_err, ChunkHeader.parseReader(elr.reader()));

    fbs = std.io.fixedBufferStream(comptime &std.mem.toBytes(std.mem.nativeToBig(u32, std.math.maxInt(u31) + 1)));
    elr = undefined;
    try std.testing.expectEqual(ChunkHeader.ParseReaderResultTag.invalid_length, ChunkHeader.parseReader(fbs.reader()));

    fbs = std.io.fixedBufferStream(comptime &std.mem.toBytes(std.mem.nativeToBig(u32, std.math.maxInt(u31))));
    elr = util.io.errorLimitedReader(fbs.reader(), 4);
    try std.testing.expectEqual(ChunkHeader.ParseReaderResultTag.no_type_err, ChunkHeader.parseReader(elr.reader()));

    fbs = std.io.fixedBufferStream(comptime &std.mem.toBytes(std.mem.nativeToBig(u32, std.math.maxInt(u31))));
    elr = undefined;
    try std.testing.expectEqual(ChunkHeader.ParseReaderResultTag.no_type_eos, ChunkHeader.parseReader(fbs.reader()));

    fbs = std.io.fixedBufferStream(comptime &std.mem.toBytes(std.mem.nativeToBig(u32, std.math.maxInt(u31))) ++ ChunkType.string(.IHDR));
    elr = undefined;
    try std.testing.expectEqual(ChunkHeader.ParseReaderResultTag.ok, ChunkHeader.parseReader(fbs.reader()));
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
        @"0" = 0,
        /// Adam7 interlace
        @"1" = 1,
        _,
    };
};
