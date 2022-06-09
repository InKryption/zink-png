const std = @import("std");
const builtin = @import("builtin");
const raw_stream = @import("raw_stream.zig");

pub const RawChunk = raw_stream.RawChunk;
pub const RawChunkStream = raw_stream.RawChunkStream;
pub const rawChunkStream = raw_stream.rawChunkStream;

pub const signature: [8]u8 = .{ 137, 80, 78, 71, 13, 10, 26, 10 };

pub const ChunkHeader = struct {
    length: u32,
    type: ChunkType,
};

pub fn chunkType(bytes: *const [4]u8) ChunkType {
    return @intToEnum(ChunkType, std.mem.readIntNative(u32, bytes));
}
pub const ChunkType = enum(u32) {
    pub const Tag = @typeInfo(ChunkType).Enum.tag_type;
    // Critical Chunk Types
    IHDR = std.mem.readIntNative(u32, "IHDR"),
    PLTE = std.mem.readIntNative(u32, "PLTE"),
    IDAT = std.mem.readIntNative(u32, "IDAT"),
    IEND = std.mem.readIntNative(u32, "IEND"),

    // Ancillary Chunk Types
    bKGD = std.mem.readIntNative(u32, "bKGD"),
    cHRM = std.mem.readIntNative(u32, "cHRM"),
    gAMA = std.mem.readIntNative(u32, "gAMA"),
    hIST = std.mem.readIntNative(u32, "hIST"),
    pHYs = std.mem.readIntNative(u32, "pHYs"),
    sBIT = std.mem.readIntNative(u32, "sBIT"),
    tEXt = std.mem.readIntNative(u32, "tEXt"),
    tIME = std.mem.readIntNative(u32, "tIME"),
    tRNS = std.mem.readIntNative(u32, "tRNS"),
    zTXt = std.mem.readIntNative(u32, "zTXt"),

    // Other Chunk Types
    _,

    pub fn intNative(self: ChunkType) Tag {
        return @enumToInt(self);
    }

    pub fn intBig(self: ChunkType) Tag {
        return std.mem.nativeToBig(Tag, self.intNative());
    }

    pub fn intLittle(self: ChunkType) Tag {
        return std.mem.nativeToLittle(Tag, self.intNative());
    }

    pub fn isCritical(self: ChunkType) bool {
        return switch (self) {
            .IHDR,
            .PLTE,
            .IDAT,
            .IEND,
            => true,

            .bKGD,
            .cHRM,
            .gAMA,
            .hIST,
            .pHYs,
            .sBIT,
            .tEXt,
            .tIME,
            .tRNS,
            .zTXt,
            => false,

            _ => false,
        };
    }

    pub fn isAncillary(self: ChunkType) bool {
        return switch (self) {
            .IHDR,
            .PLTE,
            .IDAT,
            .IEND,
            => true,

            .bKGD,
            .cHRM,
            .gAMA,
            .hIST,
            .pHYs,
            .sBIT,
            .tEXt,
            .tIME,
            .tRNS,
            .zTXt,
            => true,

            _ => false,
        };
    }
};

pub fn rawChunkBufferStream(buffer: []const u8) RawChunkBufferStream {
    return RawChunkBufferStream.init(buffer);
}

/// Wrapper over `RawChunkStream` which uses a buffer, allowing it to take advantage of the fact that it
/// can rely on the read memory remaining valid even after subsequent read calls to the reader,
/// avoiding allocations in favor of direct references to aforementioned buffer.
pub const RawChunkBufferStream = struct {
    raw_stream: Rcs,
    fbs: Fbs,

    const Rcs = RawChunkStream(Fbs.Reader);
    const Fbs = std.io.FixedBufferStream([]const u8);

    pub fn init(buffer: []const u8) RawChunkBufferStream {
        return .{
            .raw_stream = Rcs.init(undefined),
            .fbs = Fbs{
                .buffer = buffer,
                .pos = 0,
            },
        };
    }

    pub fn start(self: *RawChunkBufferStream) Rcs.StartResult {
        self.raw_stream.reader = self.fbs.reader();
        return self.raw_stream.start();
    }

    pub const NextResult = union(enum) {
        ok: RawChunk,
        no_length_bytes: Rcs.NextResult.NoLengthBytes,
        no_type_bytes: Rcs.NextResult.NoTypeBytes,
        partial_data_bytes: Rcs.NextResult.PartialDataBytes,
        no_crc_bytes: Rcs.NextResult.NoCrcBytes,
    };

    pub fn next(self: *RawChunkBufferStream) ?NextResult {
        return self.nextWithSkipSize(512);
    }

    pub fn nextWithSkipSize(self: *RawChunkBufferStream, comptime skip_buffer_bytes: usize) ?NextResult {
        var skip_buffer: [skip_buffer_bytes]u8 = undefined;
        return self.nextWithSkipBuffer(&skip_buffer);
    }

    /// Same as `RawChunkStream`, but all returned data pointers refer to the supplied buffer.
    /// Skip buffer will be used as scratch space to skip over data;
    pub fn nextWithSkipBuffer(self: *RawChunkBufferStream, skip_buffer: []u8) ?NextResult {
        const start_pos = self.fbs.pos;
        const result: Rcs.NextResult = self.raw_stream.next(.{ .skip = skip_buffer }) orelse return null;
        const end_pos = self.fbs.pos;

        const chunk_segment = self.fbs.buffer[start_pos..end_pos];

        return switch (result) {
            .ok => |info| blk: {
                const data_segment = chunk_segment[@sizeOf([2]u32) .. chunk_segment.len - @sizeOf(u32)];
                std.debug.assert(info.header.length == data_segment.len);

                break :blk NextResult{
                    .ok = RawChunk{
                        .header = info.header,
                        .p_data = data_segment.ptr,
                        .crc = info.crc,
                    },
                };
            },
            .no_length_bytes => |info| NextResult{ .no_length_bytes = info },
            .no_type_bytes => |info| NextResult{ .no_type_bytes = info },
            .out_of_mem_for_data => unreachable,
            .partial_data_bytes_no_capture => |info| blk: {
                const partial_data_segment = chunk_segment[@sizeOf([2]u32)..];
                std.debug.assert(info.bytes_len == partial_data_segment.len);

                break :blk NextResult{ .partial_data_bytes = Rcs.NextResult.PartialDataBytes{
                    .header = info.header,
                    .bytes = partial_data_segment,
                    .err = info.err,
                } };
            },
            .partial_data_bytes => unreachable,
            .no_crc_bytes_no_capture => |info| blk: {
                const data_segment = chunk_segment[@sizeOf([2]u32)..];
                std.debug.assert(info.header.length == data_segment.len);

                break :blk NextResult{ .no_crc_bytes = Rcs.NextResult.NoCrcBytes{
                    .header = info.header,
                    .p_data = data_segment.ptr,
                    .err = info.err,
                } };
            },
            .no_crc_bytes => unreachable,
        };
    }
};

test {
    std.debug.print("\n", .{});
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
                std.mem.toBytes(std.mem.nativeToBig(u8, 8)) ++ // bit depth
                std.mem.toBytes(std.mem.nativeToBig(u8, 0)) ++ // color type
                std.mem.toBytes(std.mem.nativeToBig(u8, 0)) ++ // compression method
                std.mem.toBytes(std.mem.nativeToBig(u8, 0)) ++ // filter method
                std.mem.toBytes(std.mem.nativeToBig(u8, 0)) ++ // interlace method
                // data end
                std.mem.toBytes(std.mem.nativeToBig(u32, 0x57DD52F8)), // crc
            // IDAT
            &std.mem.toBytes(std.mem.nativeToBig(u32, 17)) ++ // length
                std.mem.toBytes(ChunkType.intBig(.IDAT)) ++ // type
                [_]u8{ 8, 29, 1, 6, 0, 249, 255, 0, 255, 0, 0, 0, 255, 6, 0, 1, 255 } ++ // data
                std.mem.toBytes(std.mem.nativeToBig(u32, 0x68B6702C)), // crc
            // IEND
            &std.mem.toBytes(std.mem.nativeToBig(u32, 0)) ++ // length
                std.mem.toBytes(ChunkType.intBig(.IEND)) ++ // type
                [_]u8{} ++ // data
                std.mem.toBytes(std.mem.nativeToBig(u32, 0x68B6702C)), // crc
        }) |bytes| {
            data = data ++ bytes;
        }

        break :data data;
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    var chunk_stream = rawChunkBufferStream(data);

    try chunk_stream.start().unwrap();
    while (chunk_stream.next()) |maybe_chunk| {
        const chunk: RawChunk = switch (maybe_chunk) {
            .ok => |ok| ok,

            .no_length_bytes => |info| return info.err,
            .no_type_bytes => |info| return info.err orelse @panic("no_type_bytes"),
            .partial_data_bytes => |info| return info.err orelse @panic("partial_data_bytes"),
            .no_crc_bytes => |info| return info.err orelse @panic("no_crc_bytes"),
        };

        std.debug.print(
            \\type = '{s}'
            \\crc = 0x{X}
            \\data = [{d}]u8{{
        , .{
            @tagName(chunk.header.type),
            chunk.crc,
            chunk.header.length,
        });
        const cols = 24;
        for (chunk.data()) |byte, i| {
            if (i % cols == 0) {
                std.debug.print("\n    ", .{});
            }
            std.debug.print("{d}, ", .{byte});
        }
        if (chunk.header.length != 0) {
            std.debug.print("\n", .{});
        }
        std.debug.print("}}\n", .{});

        if (chunk.header.type == .IDAT and chunk.header.length != 0) {
            var compressed_content_stream = std.io.fixedBufferStream(chunk.data());
            var decompress_stream = try std.compress.zlib.zlibStream(arena.allocator(), compressed_content_stream.reader());
            defer decompress_stream.deinit();

            const decompressed_data: []const u8 = try decompress_stream.reader().readAllAlloc(arena.allocator(), 8);
            defer arena.allocator().free(decompressed_data);

            std.debug.print("data(decompressed) = [{d}]u8{any}\n", .{ decompressed_data.len, decompressed_data });
        }

        std.debug.print("\n", .{});
    }
}
