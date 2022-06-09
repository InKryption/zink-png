const std = @import("std");
const builtin = @import("builtin");

pub const png_signature: [8]u8 = .{ 137, 80, 78, 71, 13, 10, 26, 10 };

pub const PngChunkType = enum(u32) {
    pub const Tag = @typeInfo(PngChunkType).Enum.tag_type;
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

    _,

    pub fn intNative(self: PngChunkType) Tag {
        return @enumToInt(self);
    }

    pub fn intBig(self: PngChunkType) Tag {
        return std.mem.nativeToBig(Tag, self.intNative());
    }

    pub fn intLittle(self: PngChunkType) Tag {
        return std.mem.nativeToLittle(Tag, self.intNative());
    }

    pub fn isCritical(self: PngChunkType) bool {
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

    pub fn isAncillary(self: PngChunkType) bool {
        return switch (self) {
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
            else => false,
        };
    }
};

pub const PngChunkHeader = struct {
    length: u32,
    type: PngChunkType,
};

pub const PngRawChunk = struct {
    header: PngChunkHeader,
    p_data: [*]const u8,
    crc: u32,

    pub fn data(self: PngRawChunk) []const u8 {
        return self.p_data[0..self.header.length];
    }

    pub fn deinit(self: PngRawChunk, allocator: std.mem.Allocator) void {
        allocator.free(self.data());
    }
};

pub fn pngRawChunkStream(reader: anytype) PngRawChunkStream(@TypeOf(reader)) {
    return PngRawChunkStream(@TypeOf(reader)).init(reader);
}

pub fn PngRawChunkStream(comptime ReaderType: type) type {
    return struct {
        const Self = @This();
        reader: ReaderType,
        state: PngRawChunkStreamState,

        pub const State = PngRawChunkStreamState;
        pub const StartResult = PngRawChunkStreamStartResult(ReaderType.Error);
        pub const NextResult = PngRawChunkStreamNextResult(ReaderType.Error);

        pub fn init(reader: ReaderType) Self {
            return Self{
                .reader = reader,
                .state = .start,
            };
        }

        /// Must be called once before calling `next`.
        pub fn start(self: *Self) StartResult {
            switch (self.state) {
                .start => {
                    self.state = .end;

                    var bytes: [png_signature.len]u8 = undefined;
                    const bytes_read = self.reader.readAll(&bytes) catch |err| {
                        return StartResult{ .initial_read_fail = .{ .err = err } };
                    };
                    if (bytes_read != bytes.len) {
                        return StartResult{ .initial_read_fail = .{ .err = null } };
                    }

                    if (!std.mem.eql(u8, &bytes, &png_signature)) {
                        return .bad_png_signature;
                    }

                    self.state = .in_progress;
                },
                .in_progress => unreachable,
                .end => unreachable,
            }
            return .ok;
        }

        /// Returns 'null' when the stream ends after stepping through only whole chunks,
        /// or after failing to read a whole chunk.
        /// If no allocator is supplied, on success, the returned chunk will have an undefined pointer,
        /// and thus should not be used.
        /// Must have called `start` beforehand.
        pub fn next(self: *Self, maybe_allocator: union(enum) {
            allocator: std.mem.Allocator,
            skip: []u8,
        }) ?NextResult {
            switch (self.state) {
                .start => unreachable,
                .in_progress => {},
                .end => return null,
            }

            // Set state to 'end' by default; if nothing impedes control flow from reaching the end of the function,
            // it will be set back to `in_progress`.
            self.state = .end;

            const header: PngChunkHeader = header: {
                const length = if (util.io.readIntBigOrNull(self.reader, u32)) |maybe_int| maybe_int orelse {
                    return null;
                } else |err| return NextResult{ .no_length_bytes = NextResult.NoLengthBytes{ .err = err } };

                const @"type" = if (util.io.readIntBigOrNull(self.reader, u32)) |maybe_int| maybe_int orelse {
                    return NextResult{ .no_type_bytes = NextResult.NoTypeBytes{
                        .length = length,
                        .err = null,
                    } };
                } else |err| return NextResult{ .no_type_bytes = NextResult.NoTypeBytes{
                    .length = length,
                    .err = err,
                } };

                break :header PngChunkHeader{
                    .length = length,
                    .type = @intToEnum(PngChunkType, @"type"),
                };
            };

            const data: []const u8 = switch (maybe_allocator) {
                .allocator => |allocator| data: {
                    const data: []u8 = allocator.alloc(u8, header.length) catch |err| switch (err) {
                        error.OutOfMemory => return NextResult{ .out_of_mem_for_data = NextResult.OutOfMemForData{ .header = header } },
                    };
                    const bytes_read = self.reader.readAll(data) catch |err| {
                        allocator.free(data);
                        return NextResult{ .no_data_bytes = NextResult.NoDataBytes{
                            .header = header,
                            .err = err,
                        } };
                    };
                    std.debug.assert(bytes_read <= data.len);
                    if (bytes_read < data.len) {
                        return NextResult{ .partial_data_bytes = NextResult.PartialDataBytes{
                            .header = header,
                            .bytes = allocator.shrink(data, bytes_read),
                            .err = null,
                        } };
                    }

                    break :data data;
                },
                .skip => |buffer| data: {
                    var remaining = header.length;
                    while (remaining > 0) {
                        const amt = std.math.min(remaining, buffer.len);
                        const bytes_skipped = self.reader.readAll(buffer) catch |err| return NextResult{
                            .partial_data_bytes_no_capture = NextResult.PartialDataBytesNoCapture{
                                .header = header,
                                .bytes_len = header.length - remaining,
                                .err = err,
                            },
                        };
                        if (bytes_skipped < amt) return NextResult{
                            .partial_data_bytes_no_capture = NextResult.PartialDataBytesNoCapture{
                                .header = header,
                                .bytes_len = header.length - remaining,
                                .err = null,
                            },
                        };
                        remaining -= amt;
                    }

                    break :data undefined;
                },
            };

            const crc = if (util.io.readIntBigOrNull(self.reader, u32)) |maybe_int| maybe_int orelse return switch (std.meta.activeTag(maybe_allocator)) {
                .skip => NextResult{ .no_crc_bytes_no_capture = NextResult.NoCrcBytesNoCapture{
                    .header = header,
                    .err = null,
                } },
                .allocator => NextResult{ .no_crc_bytes = NextResult.NoCrcBytes{
                    .header = header,
                    .p_data = data.ptr,
                    .err = null,
                } },
            } else |err| return switch (std.meta.activeTag(maybe_allocator)) {
                .skip => NextResult{ .no_crc_bytes_no_capture = NextResult.NoCrcBytesNoCapture{
                    .header = header,
                    .err = err,
                } },
                .allocator => NextResult{ .no_crc_bytes = NextResult.NoCrcBytes{
                    .header = header,
                    .p_data = data.ptr,
                    .err = err,
                } },
            };

            // Having reached the end of the function without any issues, set state back to `in_progress`.
            self.state = .in_progress;

            return NextResult{
                .ok = PngRawChunk{
                    .header = header,
                    .p_data = data.ptr,
                    .crc = crc,
                },
            };
        }
    };
}

const PngRawChunkStreamState = enum {
    start,
    in_progress,
    end,
};
const PngRawChunkStreamStartResultTag = enum {
    ok,
    initial_read_fail,
    bad_png_signature,
};
const PngRawChunkStreamNextResultTag = enum {
    ok,
    no_length_bytes,
    no_type_bytes,
    no_data_bytes,
    out_of_mem_for_data,
    partial_data_bytes_no_capture,
    partial_data_bytes,
    no_crc_bytes_no_capture,
    no_crc_bytes,
};

fn PngRawChunkStreamStartResult(comptime ReadError: type) type {
    return union(PngRawChunkStreamStartResultTag) {
        const Self = @This();
        ok,
        initial_read_fail: InitialReadFail,
        bad_png_signature,

        pub const InitialReadFail = struct {
            /// 'null' if stream ended.
            err: ?ReadError,
        };

        pub const Error = ReadError || error{
            EndOfStream,
            BadPngSignature,
        };

        pub fn unwrap(self: Self) Error!void {
            return switch (self) {
                .ok => {},
                .initial_read_fail => |info| info.err orelse error.EndOfStream,
                .bad_png_signature => error.BadPngSignature,
            };
        }
    };
}
fn PngRawChunkStreamNextResult(comptime ReadError: type) type {
    return union(PngRawChunkStreamNextResultTag) {
        const Self = @This();
        ok: PngRawChunk,
        no_length_bytes: NoLengthBytes,
        no_type_bytes: NoTypeBytes,
        no_data_bytes: NoDataBytes,
        out_of_mem_for_data: OutOfMemForData,
        partial_data_bytes_no_capture: PartialDataBytesNoCapture,
        partial_data_bytes: PartialDataBytes,
        no_crc_bytes_no_capture: NoCrcBytesNoCapture,
        no_crc_bytes: NoCrcBytes,

        pub const Error = ReadError;

        pub const NoLengthBytes = struct {
            err: Error,
        };
        pub const NoTypeBytes = struct {
            length: u32,
            /// 'null' if stream ended.
            err: ?Error,
        };
        pub const NoDataBytes = HeaderMaybeErrOnly;
        pub const OutOfMemForData = struct {
            header: PngChunkHeader,
        };
        pub const PartialDataBytesNoCapture = struct {
            header: PngChunkHeader,
            bytes_len: usize,
            /// 'null' if stream ended.
            err: ?Error,
        };
        pub const PartialDataBytes = struct {
            header: PngChunkHeader,
            bytes: []const u8,
            /// 'null' if stream ended.
            err: ?Error,
        };
        pub const NoCrcBytes = struct {
            header: PngChunkHeader,
            p_data: [*]const u8,
            /// 'null' if stream ended.
            err: ?Error,
        };
        pub const NoCrcBytesNoCapture = HeaderMaybeErrOnly;

        const HeaderMaybeErrOnly = struct {
            header: PngChunkHeader,
            /// 'null' if stream ended.
            err: ?Error,
        };
    };
}

test {
    std.debug.print("\n", .{});
    const data: []const u8 = comptime data: {
        var data: []const u8 = "";

        for ([_][]const u8{
            // signature
            &png_signature,
            // IHDR
            &std.mem.toBytes(std.mem.nativeToBig(u32, 13)) ++ // length
                std.mem.toBytes(PngChunkType.intBig(.IHDR)) ++ // type
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
                std.mem.toBytes(PngChunkType.intBig(.IDAT)) ++ // type
                [_]u8{ 8, 29, 1, 6, 0, 249, 255, 0, 255, 0, 0, 0, 255, 6, 0, 1, 255 } ++ // data
                std.mem.toBytes(std.mem.nativeToBig(u32, 0x68B6702C)), // crc
            // IEND
            &std.mem.toBytes(std.mem.nativeToBig(u32, 0)) ++ // length
                std.mem.toBytes(PngChunkType.intBig(.IEND)) ++ // type
                [_]u8{} ++ // data
                std.mem.toBytes(std.mem.nativeToBig(u32, 0x68B6702C)), // crc
        }) |bytes| {
            data = data ++ bytes;
        }

        break :data data;
    };
    var data_stream = std.io.fixedBufferStream(data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    var chunk_stream = pngRawChunkStream(data_stream.reader());

    try chunk_stream.start().unwrap();
    while (chunk_stream.next(.{ .allocator = arena.allocator() })) |maybe_chunk| {
        const chunk: PngRawChunk = switch (maybe_chunk) {
            .ok => |ok| ok,

            .no_length_bytes => |info| return info.err,
            .no_type_bytes => |info| return info.err orelse @panic("no_type_bytes"),
            .no_data_bytes => |info| return info.err orelse @panic("no_data_bytes"),
            .out_of_mem_for_data => return error.OutOfMemory,
            .partial_data_bytes_no_capture => |info| return info.err orelse @panic("partial_data_bytes_no_capture"),
            .partial_data_bytes => |info| return info.err orelse @panic("partial_data_bytes"),
            .no_crc_bytes_no_capture => |info| return info.err orelse @panic("no_crc_bytes_no_capture"),
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

const util = struct {
    const io = struct {
        pub fn readBytesNoEof(reader: anytype, comptime num_bytes: usize) @TypeOf(reader).Error!?[num_bytes]u8 {
            var buffer: [num_bytes]u8 = undefined;
            const bytes_read = try reader.readAll(&buffer);
            return if (bytes_read < buffer.len) null else buffer;
        }

        /// Tries to read a BE integer from the stream; if the stream ends before supplying enough bytes
        /// for such an integer, returns null.
        pub fn readIntBigOrNull(reader: anytype, comptime T: type) @TypeOf(reader).Error!?T {
            const byte_count = (@typeInfo(T).Int.bits + 7) / 8;
            const bytes = (try readBytesNoEof(reader, byte_count)) orelse return null;
            return std.mem.readIntBig(T, &bytes);
        }
    };
};
