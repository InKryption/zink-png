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

    pub fn str(self: ChunkType) [4]u8 {
        return std.mem.toBytes(self.int());
    }

    pub fn isValid(self: ChunkType) bool {
        for (self.str()) |byte| {
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
        return (self.str()[byte_index] & 32) != 0;
    }
};

pub const ChunkHeader = struct {
    length: u31,
    type: ChunkType,

    pub fn parseBytes(bytes: *const [@sizeOf(u32) * 2]u8) ChunkHeader {
        return switch (ChunkHeader.parseBuffer(bytes)) {
            .ok => |value| value,
            else => unreachable,
        };
    }

    pub const ParseBufferResult = union(enum) {
        ok: ChunkHeader,
        no_type: NoType,
        no_length: NoLength,

        pub const NoType = struct { length: u32 };
        pub const NoLength = struct {};
    };

    pub fn parseBuffer(buffer: []const u8) ParseBufferResult {
        var fbs = std.io.fixedBufferStream(buffer);
        return switch (ChunkHeader.parseReader(fbs.reader())) {
            .ok => |value| ParseBufferResult{ .ok = value },
            .no_type_eos => |info| ParseBufferResult{ .no_type = .{ .length = info.length } },
            .no_type_err => unreachable,
            .no_length_eos => ParseBufferResult{ .no_length = .{} },
            .no_length_err => unreachable,
        };
    }

    pub fn ParseReaderResult(comptime ReaderError: type) type {
        return union(enum) {
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
            pub const NoTypeEos = struct { length: u31 };
            pub const NoTypeErr = struct { length: u31, err: ReadError };
            pub const InvalidLength = struct { length: u32 };
            pub const NoLengthEos = struct {};
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
