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

pub fn rawChunkStream(reader: anytype) RawChunkStream(@TypeOf(reader)) {
    return RawChunkStream(@TypeOf(reader)).init(reader);
}
pub fn RawChunkStream(comptime Reader: type) type {
    return struct {
        const Self = @This();
        reader: Reader,
        state: State,

        pub fn init(reader: Reader) Self {
            return Self{
                .reader = reader,
                .state = .begin,
            };
        }

        pub const StartError = error{ NoBytes, IncompleteSignature, InvalidSignature };
        pub fn start(self: *Self) Reader.Error!StartError!void {
            switch (self.state) {
                .begin => {},
                .awaiting_header => unreachable,
                .awaiting_data => unreachable,
                .end => unreachable,
            }

            const Result = StartError!void;
            errdefer self.state = .end;

            const start_bytes = try util.io.readBoundedArray(self.reader, signature.len);

            if (start_bytes.len == 0) return util.as(Result, error.NoBytes);
            if (start_bytes.len < signature.len) return util.as(Result, error.IncompleteSignature);

            std.debug.assert(start_bytes.len == signature.len);
            if (!std.mem.eql(u8, start_bytes.constSlice(), signature[0..])) {
                return util.as(Result, error.InvalidSignature);
            }

            self.state = .awaiting_header;
        }

        const GetHeaderResult = union(enum) {
            /// success
            ok: ChunkHeader,
            /// encountered end of stream while trying to read type
            no_type_eos: NoTypeEos,
            /// encountered error while trying to read type
            no_type_err: NoTypeErr,
            /// encountered an invalid length
            invalid_length: InvalidLength,
            /// encountered error while trying to read length
            no_length_err: NoLengthErr,

            const ReadError = util.NormalizedErrorSet(Reader.Error);
            pub const NoTypeEos = ChunkHeader.ParseReaderResult(ReadError).NoTypeEos;
            pub const NoTypeErr = ChunkHeader.ParseReaderResult(ReadError).NoTypeErr;
            pub const InvalidLength = ChunkHeader.ParseReaderResult(ReadError).InvalidLength;
            pub const NoLengthErr = ChunkHeader.ParseReaderResult(ReadError).NoLengthErr;
        };
        pub fn getHeader(self: *Self) ?GetHeaderResult {
            switch (self.state) {
                .begin => unreachable,
                .awaiting_header => {},
                .awaiting_data => unreachable,
                .end => return null,
            }

            const result = ChunkHeader.parseReader(self.reader);
            self.state = switch (result) {
                .ok => |header| State{ .awaiting_data = .{ .header = header } },
                .no_type_eos,
                .no_type_err,
                .invalid_length,
                .no_length_eos,
                .no_length_err,
                => .end,
            };
            return switch (result) {
                .ok => |header| GetHeaderResult{ .ok = header },
                .no_type_eos => |info| GetHeaderResult{ .no_type_eos = info },
                .no_type_err => |info| GetHeaderResult{ .no_type_err = info },
                .invalid_length => |info| GetHeaderResult{ .invalid_length = info },
                .no_length_eos => null,
                .no_length_err => |info| GetHeaderResult{ .no_length_err = info },
            };
        }

        pub fn getDataAndCrcWithBuffer(self: *Self, writer: anytype, intermediate_buffer: []u8) Reader.Error!@TypeOf(writer).Error!?u32 {
            const Result = @TypeOf(writer).Error!?u32;

            std.debug.assert(intermediate_buffer.len >= 1);
            const ctx: State.AwaitingData = switch (self.state) {
                .begin => unreachable,
                .awaiting_header => unreachable,
                .awaiting_data => |ctx| ctx,
                .end => unreachable,
            };
            errdefer self.state = .end;

            var i: u32 = 0;
            while (i != ctx.header.length) {
                const buf_len = std.math.min(intermediate_buffer.len, ctx.header.length - i);
                const amt = try self.reader.read(intermediate_buffer[0..buf_len]);

                writer.writeAll(intermediate_buffer[0..amt]) catch |err| {
                    self.state = .end;
                    return util.as(Result, err);
                };

                i += amt;
                if (amt == 0) {
                    self.state = .end;
                    return null;
                }
            }

            return util.io.readIntBigOrNull(self.reader, u32) catch |err| return err;
        }

        const State = union(enum) {
            begin,
            awaiting_header,
            awaiting_data: AwaitingData,
            end,

            const AwaitingData = struct { header: ChunkHeaderwell };
        };
    };
}

test "RawChunkStream.start + FixedBufferStream" {
    var fbs: std.io.FixedBufferStream([]const u8) = undefined;
    var rcs: RawChunkStream(@TypeOf(fbs).Reader) = undefined;

    fbs = std.io.fixedBufferStream(&[_]u8{});
    rcs = rawChunkStream(fbs.reader());
    try std.testing.expectError(error.NoBytes, rcs.start() catch |err| switch (err) {});

    fbs = std.io.fixedBufferStream(signature[0 .. signature.len - 1]);
    rcs = rawChunkStream(fbs.reader());
    try std.testing.expectError(error.IncompleteSignature, rcs.start() catch |err| switch (err) {});

    fbs = std.io.fixedBufferStream(signature[0 .. signature.len - 1] ++ [_]u8{signature[signature.len - 1] +% 1});
    rcs = rawChunkStream(fbs.reader());
    try std.testing.expectError(error.InvalidSignature, rcs.start() catch |err| switch (err) {});
}

test "RawChunkStream.start + fallible reader" {
    var fbs: std.io.FixedBufferStream([]const u8) = undefined;
    var elr: util.io.ErrorLimitedReader(@TypeOf(fbs).Reader) = undefined;
    var rcs: RawChunkStream(@TypeOf(elr).Reader) = undefined;

    fbs = std.io.fixedBufferStream(signature[0..]);
    elr = util.io.errorLimitedReader(fbs.reader(), signature.len - 1);
    rcs = rawChunkStream(elr.reader());
    try std.testing.expectError(error.AttemptedToReadMoreThanLimit, rcs.start());
}

test "RawChunkStream.getHeader" {
    var fbs: std.io.FixedBufferStream([]const u8) = undefined;
    var rcs: RawChunkStream(@TypeOf(fbs).Reader) = undefined;

    fbs = std.io.fixedBufferStream(signature[0..]);
    rcs = rawChunkStream(fbs.reader());
    rcs.start() catch |err| switch (err) {} catch @panic("That wasn't supposed to happen");
    try std.testing.expectEqual(@as(?@TypeOf(rcs).GetHeaderResult, null), rcs.getHeader());

    fbs = std.io.fixedBufferStream(signature[0..] ++
        std.mem.toBytes(std.mem.nativeToBig(u32, 13)) ++
        std.mem.toBytes(ChunkType.intBig(.IHDR)));
    rcs = rawChunkStream(fbs.reader());
    rcs.start() catch |err| switch (err) {} catch @panic("That wasn't supposed to happen");
    try std.testing.expectEqual(@as(?@TypeOf(rcs).GetHeaderResult, @TypeOf(rcs).GetHeaderResult{ .ok = ChunkHeader{ .length = 13, .type = .IHDR } }), rcs.getHeader());
}
