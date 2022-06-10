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
    length: u32,
    type: ChunkType,

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

    pub const ParseBufferResult = union(enum) {
        ok: ChunkHeader,
        no_type: NoType,
        no_length: NoLength,

        pub const NoType = struct { length: u32 };
        pub const NoLength = struct {};
    };

    /// Returns null if the stream ends before returning the required number of bytes for a chunk header.
    pub fn parseReader(reader: anytype) ParseReaderResult(util.MemoizedErrorSet(@TypeOf(reader).Error)) {
        const PResult = ParseReaderResult(util.MemoizedErrorSet(@TypeOf(reader).Error));

        const length = util.io.readIntBigOrNull(reader, u32) catch |err| {
            return PResult{ .no_length_err = .{ .err = err } };
        } orelse return PResult{ .no_length_eos = .{} };

        const type_value = util.io.readIntBigOrNull(reader, u32) catch |err| {
            return PResult{ .no_type_err = .{ .err = err, .length = length } };
        } orelse return PResult{ .no_type_eos = .{ .length = length } };

        return PResult{ .ok = ChunkHeader{
            .length = length,
            .type = @intToEnum(ChunkType, type_value),
        } };
    }

    const ParseReaderResultTag = enum {
        /// success
        ok,
        /// encountered end of stream while trying to read type
        no_type_eos,
        /// encountered error while trying to read type
        no_type_err,
        /// encountered end of stream while trying to read length
        no_length_eos,
        /// encountered error while trying to read length
        no_length_err,
    };
    pub fn ParseReaderResult(comptime ReaderError: type) type {
        return union(ParseReaderResultTag) {
            const Self = @This();
            ok: ChunkHeader,
            no_length_err: NoLengthErr,
            no_length_eos: NoLengthEos,
            no_type_err: NoTypeErr,
            no_type_eos: NoTypeEos,

            pub const ReadError = ReaderError;
            pub const NoLengthErr = struct { err: ReadError };
            pub const NoLengthEos = struct {};
            pub const NoTypeErr = struct { length: u32, err: ReadError };
            pub const NoTypeEos = struct { length: u32 };
        };
    }
};

pub const RawChunk = struct {
    header: ChunkHeader,
    p_data: [*]const u8,
    crc: u32,

    pub fn data(self: RawChunk) []const u8 {
        return self.p_data[0..self.header.length];
    }

    pub fn deinit(self: RawChunk, allocator: std.mem.Allocator) void {
        allocator.free(self.data());
    }
};

pub const RawChunkStreamMemory = struct {
    src: []const u8,
    state: State,
    index: usize,

    pub fn init(buffer: []const u8) RawChunkStreamMemory {
        return RawChunkStreamMemory{
            .src = buffer,
            .state = .begin,
            .index = 0,
        };
    }

    pub const StartError = error{ NoPngSignature, BadPngSignature };
    pub fn start(self: *RawChunkStreamMemory) StartError!void {
        switch (self.state) {
            .begin => {
                errdefer self.state = .end;

                std.debug.assert(self.index == 0);
                self.index += signature.len;

                if (self.src.len < self.index) {
                    return error.NoPngSignature;
                }
                if (!std.mem.eql(u8, self.src[0..self.index], &signature)) {
                    return error.BadPngSignature;
                }

                self.state = .in_progress;
            },
            .in_progress => unreachable,
            .end => unreachable,
        }
    }

    pub const NextResult = union(enum) {
        ok: RawChunk,
        no_crc: NoCrc,
        partial_data: PartialData,
        no_type: NoType,
        no_length: NoLength,

        pub const UnwrapError = error{ NoCrc, PartialData, NoType, NoLength };
        pub fn unwrap(self: NextResult) UnwrapError!RawChunk {
            return switch (self) {
                .ok => |value| value,
                .no_crc => error.NoCrc,
                .partial_data => error.PartialData,
                .no_type => error.NoType,
                .no_length => error.NoLength,
            };
        }

        pub const NoCrc = struct { header: ChunkHeader, p_data: [*]const u8 };
        pub const PartialData = struct { header: ChunkHeader, partial_data: []const u8 };
        pub const NoType = struct { length: u32 };
        pub const NoLength = struct {};
    };
    pub fn next(self: *RawChunkStreamMemory) ?NextResult {
        switch (self.state) {
            .begin => unreachable,
            .in_progress => {},
            .end => return null,
        }

        self.state = .end;
        if (self.index == self.src.len) {
            return null;
        }

        const header: ChunkHeader = switch (ChunkHeader.parseBuffer(self.src[self.index..])) {
            .ok => |value| value,
            .no_type => |info| return NextResult{ .no_type = .{ .length = info.length } },
            .no_length => return NextResult{ .no_length = .{} },
        };
        self.index += @sizeOf([2]u32);

        const data: []const u8 = data: {
            const data_start = self.index;
            self.index += header.length;

            if (self.src.len < self.index) {
                return NextResult{ .partial_data = NextResult.PartialData{
                    .header = header,
                    .partial_data = self.src[data_start..],
                } };
            }

            const data = self.src[data_start..self.index];
            std.debug.assert(data.len == header.length);

            break :data data;
        };

        const crc: u32 = crc: {
            const crc_start = self.index;
            self.index += @sizeOf(u32);

            if (self.src.len < self.index) {
                return NextResult{ .no_crc = NextResult.NoCrc{
                    .header = header,
                    .p_data = data.ptr,
                } };
            }

            const crc = std.mem.readIntBig(u32, self.src[crc_start..][0..@sizeOf(u32)]);
            break :crc crc;
        };

        self.state = .in_progress;

        return NextResult{ .ok = RawChunk{
            .header = header,
            .p_data = data.ptr,
            .crc = crc,
        } };
    }

    const State = enum {
        begin,
        in_progress,
        end,
    };
};

test {
    const data: []const u8 = comptime data: {
        var data: []const u8 = "";

        for ([_][]const u8{
            // signature
            &signature,

            // IHDR
            &std.mem.toBytes(std.mem.nativeToBig(u32, 13)) ++ // length
                std.mem.toBytes(ChunkType.intBig(.IHDR)) ++ // type
                // data start
                std.mem.toBytes(std.mem.nativeToBig(u32, 2)) ++ // width
                std.mem.toBytes(std.mem.nativeToBig(u32, 2)) ++ // height
                [_]u8{
                8, // bit depth
                0, // color type
                0, // compression method
                0, // filter method
                0, // interlace method
            } ++ std.mem.toBytes(std.mem.nativeToBig(u32, 0x57DD52F8)), // crc

            // IDAT
            &std.mem.toBytes(std.mem.nativeToBig(u32, 17)) ++ // length
                std.mem.toBytes(ChunkType.intBig(.IDAT)) ++ // type
                [17]u8{ 8, 29, 1, 6, 0, 249, 255, 0, 255, 0, 0, 0, 255, 6, 0, 1, 255 } ++ // data
                std.mem.toBytes(std.mem.nativeToBig(u32, 0x68B6702C)), // crc

            // IEND
            &std.mem.toBytes(std.mem.nativeToBig(u32, 0)) ++ // length
                std.mem.toBytes(ChunkType.intBig(.IEND)) ++ // type
                [0]u8{} ++ // data
                std.mem.toBytes(std.mem.nativeToBig(u32, 0xAE426082)), // crc
        }) |bytes| {
            data = data ++ bytes;
        }

        break :data data;
    };

    var rcsm = RawChunkStreamMemory.init(data);

    try rcsm.start();
    if (rcsm.next()) |maybe_chunk| {
        const chunk = try maybe_chunk.unwrap();

        try std.testing.expect(chunk.header.type.isValid());
        try std.testing.expectEqual(ChunkType.from("IHDR"), chunk.header.type);
        try std.testing.expectEqualStrings(
            &std.mem.toBytes(std.mem.nativeToBig(u32, 2)) ++
                std.mem.toBytes(std.mem.nativeToBig(u32, 2)) ++
                [_]u8{ 8, 0, 0, 0, 0 },
            chunk.data(),
        );

        var crc_hasher = std.hash.Crc32.init();
        crc_hasher.update(&std.mem.toBytes(chunk.header.type.intBig()));
        crc_hasher.update(chunk.data());
        try std.testing.expectEqual(crc_hasher.final(), chunk.crc);
    } else return error.UnexpectedNullChunk;

    if (rcsm.next()) |maybe_chunk| {
        const chunk = try maybe_chunk.unwrap();
        try std.testing.expect(chunk.header.type.isValid());
        try std.testing.expectEqual(ChunkType.from("IDAT"), chunk.header.type);
        try std.testing.expectEqualSlices(
            u8,
            &[17]u8{ 8, 29, 1, 6, 0, 249, 255, 0, 255, 0, 0, 0, 255, 6, 0, 1, 255 },
            chunk.data(),
        );

        var crc_hasher = std.hash.Crc32.init();
        crc_hasher.update(&std.mem.toBytes(chunk.header.type.intBig()));
        crc_hasher.update(chunk.data());
        try std.testing.expectEqual(crc_hasher.final(), chunk.crc);

        var chunk_data_stream = std.io.fixedBufferStream(chunk.data());
        var zlib_stream = try std.compress.zlib.zlibStream(std.testing.allocator, chunk_data_stream.reader());
        defer zlib_stream.deinit();

        const filtered_contents = try zlib_stream.reader().readBytesNoEof(6);
        try std.testing.expectEqualSlices(
            u8,
            &[_]u8{
                0, 255, 000,
                0, 000, 255,
            },
            &filtered_contents,
        );

        try std.testing.expectEqual(@as(usize, 0), blk: {
            var skip_buff: [1]u8 = undefined;
            break :blk try zlib_stream.reader().readAll(&skip_buff);
        });
    } else return error.UnexpectedNullChunk;

    if (rcsm.next()) |maybe_chunk| {
        const chunk = try maybe_chunk.unwrap();
        try std.testing.expect(chunk.header.type.isValid());
        try std.testing.expectEqual(ChunkType.from("IEND"), chunk.header.type);
        try std.testing.expectEqualSlices(u8, &[0]u8{}, chunk.data());

        var crc_hasher = std.hash.Crc32.init();
        crc_hasher.update(&std.mem.toBytes(chunk.header.type.intBig()));
        crc_hasher.update(chunk.data());
        try std.testing.expectEqual(crc_hasher.final(), chunk.crc);
    } else return error.UnexpectedNullChunk;

    try std.testing.expectEqual(@as(?RawChunkStreamMemory.NextResult, null), rcsm.next());
}
