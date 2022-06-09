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
        no_type_bytes: NoTypeBytes,
        partial_data_bytes: PartialDataBytes,
        no_crc_bytes: NoCrcBytes,

        pub const NoTypeBytes = struct { length: u32 };
        pub const PartialDataBytes = struct { header: ChunkHeader, bytes: []const u8 };
        pub const NoCrcBytes = struct { header: ChunkHeader, p_data: [*]const u8 };

        pub const UnwrapError = error{ NoTypeBytes, IncompleteData, NoCrcBytes };
        pub fn unwrap(self: NextResult) UnwrapError!RawChunk {
            return switch (self) {
                .ok => |chunk| chunk,
                .no_type_bytes => error.NoTypeBytes,
                .partial_data_bytes => error.IncompleteData,
                .no_crc_bytes => error.NoCrcBytes,
            };
        }
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
                        .p_data = std.mem.span(data_segment).ptr,
                        .crc = info.crc,
                    },
                };
            },
            .no_length_bytes => unreachable,
            .no_type_bytes => |info| blk: {
                std.debug.assert(info.err == null);
                break :blk NextResult{ .no_type_bytes = NextResult.NoTypeBytes{
                    .length = info.length,
                } };
            },
            .out_of_mem_for_data => unreachable,
            .partial_data_bytes_no_capture => |info| blk: {
                const partial_data_segment = chunk_segment[@sizeOf([2]u32)..];
                std.debug.assert(info.bytes_len == partial_data_segment.len);
                std.debug.assert(info.err == null);

                break :blk NextResult{ .partial_data_bytes = NextResult.PartialDataBytes{
                    .header = info.header,
                    .bytes = std.mem.span(partial_data_segment),
                } };
            },
            .partial_data_bytes => unreachable,
            .no_crc_bytes_no_capture => |info| blk: {
                const data_segment = chunk_segment[@sizeOf([2]u32)..];
                std.debug.assert(info.header.length == data_segment.len);
                std.debug.assert(info.err == null);

                break :blk NextResult{ .no_crc_bytes = NextResult.NoCrcBytes{
                    .header = info.header,
                    .p_data = std.mem.span(data_segment).ptr,
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
                [17]u8{ 8, 29, 1, 6, 0, 249, 255, 0, 255, 0, 0, 0, 255, 6, 0, 1, 255 } ++ // data
                std.mem.toBytes(std.mem.nativeToBig(u32, 0x68B6702C)), // crc
            // IEND
            &std.mem.toBytes(std.mem.nativeToBig(u32, 0)) ++ // length
                std.mem.toBytes(ChunkType.intBig(.IEND)) ++ // type
                [0]u8{} ++ // data
                std.mem.toBytes(std.mem.nativeToBig(u32, 0x68B6702C)), // crc
        }) |bytes| {
            data = data ++ bytes;
        }

        break :data data;
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    comptime var chunk_stream = rawChunkBufferStream(data);

    comptime try chunk_stream.start().unwrap();
    inline while (comptime chunk_stream.next()) |maybe_chunk| {
        const chunk: RawChunk = comptime try maybe_chunk.unwrap();
        const chunk_data = comptime chunk.data();

        std.debug.print(
            \\type = '{s}'
            \\crc = 0x{X}
            \\data = [{d}]u8{any}
            \\
            \\
        , .{
            @tagName(chunk.header.type),
            chunk.crc,
            chunk_data.len,
            chunk_data,
        });

        if (comptime chunk.header.type == .IDAT and chunk.header.length != 0) {
            
            var compressed_content_stream = comptime std.io.fixedBufferStream(chunk_data);
            // The issue manifests here, although it originates further back probably;
            // `chunk_data` as passed seems to start 8 bytes before where it should, causing
            // the ZlibStream to detect the first byte as being incorrect - this doesn't happen at runtime
            // where the `chunk_data` slice is correct.
            var decompress_stream = try std.compress.zlib.zlibStream(arena.allocator(), compressed_content_stream.reader());
            defer decompress_stream.deinit();

            const decompressed_data: []const u8 = try decompress_stream.reader().readAllAlloc(arena.allocator(), 8);
            defer arena.allocator().free(decompressed_data);

            std.debug.print("data(decompressed) = [{d}]u8{any}\n", .{ decompressed_data.len, decompressed_data });
        }

        std.debug.print("\n", .{});
    }
}
