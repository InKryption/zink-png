const std = @import("std");
const builtin = @import("builtin");

const png_signature: [8]u8 = .{ 137, 80, 78, 71, 13, 10, 26, 10 };

pub const PngChunkType = enum(u32) {
    pub const Tag = @typeInfo(PngChunkType).Enum.tag_type;
    // Critical Chunk Types
    IHDR = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "IHDR")),
    PLTE = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "PLTE")),
    IDAT = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "IDAT")),
    IEND = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "IEND")),

    // Ancillary Chunk Types
    bKGD = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "bKGD")),
    cHRM = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "cHRM")),
    gAMA = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "gAMA")),
    hIST = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "hIST")),
    pHYs = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "pHYs")),
    sBIT = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "sBIT")),
    tEXt = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "tEXt")),
    tIME = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "tIME")),
    tRNS = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "tRNS")),
    zTXt = std.mem.nativeToBig(u32, std.mem.readIntNative(u32, "zTXt")),

    _,

    pub fn intNative(self: PngChunkType) Tag {
        return std.mem.bigToNative(Tag, @enumToInt(self));
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
        state: State,

        pub fn init(reader: ReaderType) Self {
            return Self{
                .reader = reader,
                .state = .start,
            };
        }

        const State = enum {
            start,
            in_progress,
            end,
        };

        pub const Result = PngRawChunkStreamResult(ReaderType.Error);
        /// Returns 'null' when the stream ends after stepping through only whole chunks,
        /// or after failing to read a whole chunk.
        /// If no allocator is supplied, on success, the returned chunk will have an undefined pointer,
        /// and thus should not be used.
        pub fn next(self: *Self, maybe_allocator: union(enum) {
            allocator: std.mem.Allocator,
            skip: []u8,
        }) ?Result {
            switch (self.state) {
                .start => {
                    self.state = .end;

                    var bytes: [png_signature.len]u8 = undefined;
                    if (util.io.readNoEof(self.reader, &bytes)) |maybe_bytes_read|
                        (if (maybe_bytes_read != null) return Result{ .initial_read_fail = .{ .err = null } })
                    else |err| {
                        return Result{ .initial_read_fail = .{ .err = err } };
                    }

                    if (!std.mem.eql(u8, &bytes, &png_signature)) {
                        return .bad_png_signature;
                    }

                    self.state = .in_progress;
                },
                .in_progress => {},
                .end => return null,
            }

            // Set state to 'end' by default; if nothing impedes control flow from reaching the end of the function,
            // it will be set back to `in_progress`.
            self.state = .end;

            const header: PngChunkHeader = header: {
                const length = if (util.io.readIntBigOrNull(self.reader, u32)) |maybe_int| maybe_int orelse {
                    return null;
                } else |err| return Result{ .no_length_bytes = Result.NoLengthBytes{ .err = err } };

                const @"type" = if (util.io.readIntBigOrNull(self.reader, u32)) |maybe_int| maybe_int orelse {
                    return Result{ .no_type_bytes = Result.NoTypeBytes{
                        .length = length,
                        .err = null,
                    } };
                } else |err| return Result{ .no_type_bytes = Result.NoTypeBytes{
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
                        error.OutOfMemory => return .out_of_memory,
                    };
                    const bytes_read = self.reader.readAll(data) catch |err| {
                        allocator.free(data);
                        return Result{ .no_data_bytes = Result.NoDataBytes{
                            .header = header,
                            .err = err,
                        } };
                    };
                    std.debug.assert(bytes_read <= data.len);
                    if (bytes_read < data.len) {
                        return Result{ .partial_data_bytes = Result.PartialDataBytes{
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
                        const bytes_skipped = self.reader.readAll(buffer) catch |err| return Result{
                            .partial_data_bytes_no_capture = Result.PartialDataBytesNoCapture{
                                .header = header,
                                .bytes_len = header.length - remaining,
                                .err = err,
                            },
                        };
                        if (bytes_skipped < amt) return Result{
                            .partial_data_bytes_no_capture = Result.PartialDataBytesNoCapture{
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
                .skip => Result{ .no_crc_bytes_no_capture = Result.NoCrcBytesNoCapture{
                    .header = header,
                    .err = null,
                } },
                .allocator => Result{ .no_crc_bytes = Result.NoCrcBytes{
                    .header = header,
                    .p_data = data.ptr,
                    .err = null,
                } },
            } else |err| return switch (std.meta.activeTag(maybe_allocator)) {
                .skip => Result{ .no_crc_bytes_no_capture = Result.NoCrcBytesNoCapture{
                    .header = header,
                    .err = err,
                } },
                .allocator => Result{ .no_crc_bytes = Result.NoCrcBytes{
                    .header = header,
                    .p_data = data.ptr,
                    .err = err,
                } },
            };

            // Having reached the end of the function without any issues, set state back to `in_progress`.
            self.state = .in_progress;

            return Result{
                .ok = PngRawChunk{
                    .header = header,
                    .p_data = data.ptr,
                    .crc = crc,
                },
            };
        }
    };
}

pub fn PngRawChunkStreamResult(comptime Error: type) type {
    return union(enum) {
        ok: PngRawChunk,

        out_of_memory,

        initial_read_fail: InitialReadFail,
        bad_png_signature,

        no_length_bytes: NoLengthBytes,
        no_type_bytes: NoTypeBytes,
        no_data_bytes: NoDataBytes,
        partial_data_bytes_no_capture: PartialDataBytesNoCapture,
        partial_data_bytes: PartialDataBytes,
        no_crc_bytes_no_capture: NoCrcBytesNoCapture,
        no_crc_bytes: NoCrcBytes,

        pub const InitialReadFail = struct {
            /// 'null' if stream ended.
            err: ?Error,
        };
        pub const NoLengthBytes = struct {
            err: Error,
        };
        pub const NoTypeBytes = struct {
            length: u32,
            /// 'null' if stream ended.
            err: ?Error,
        };
        pub const NoDataBytes = struct {
            header: PngChunkHeader,
            /// 'null' if stream ended.
            err: ?Error,
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
        pub const NoCrcBytesNoCapture = struct {
            header: PngChunkHeader,
            /// 'null' if stream ended.
            err: ?Error,
        };
    };
}

test {
    std.debug.print("\n", .{});
    const data: []const u8 = @embedFile("img1.png");
    var data_stream = std.io.fixedBufferStream(data);

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    var chunk_stream = pngRawChunkStream(data_stream.reader());
    while (chunk_stream.next(.{ .allocator = arena.allocator() })) |maybe_chunk| {
        const chunk: PngRawChunk = switch (maybe_chunk) {
            .ok => |ok| ok,

            .out_of_memory => return error.OutOfMemory,

            .initial_read_fail => |info| return info.err orelse @panic("initial_read_fail"),
            .bad_png_signature => @panic("bad_png_signature"),

            .no_length_bytes => |info| return info.err,
            .no_type_bytes => |info| return info.err orelse @panic("no_type_bytes"),
            .no_data_bytes => |info| return info.err orelse @panic("no_data_bytes"),
            .partial_data_bytes_no_capture => |info| return info.err orelse @panic("partial_data_bytes_no_capture"),
            .partial_data_bytes => |info| return info.err orelse @panic("partial_data_bytes"),
            .no_crc_bytes_no_capture => |info| return info.err orelse @panic("no_crc_bytes_no_capture"),
            .no_crc_bytes => |info| return info.err orelse @panic("no_crc_bytes"),
        };

        std.debug.print(
            \\type = '{s}'
            \\crc = {x}
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
        std.debug.print("}}\n\n", .{});
    }
}

const util = struct {
    const io = struct {
        /// Attempts to read bytes into the entire buffer; if the stream ends before filling the provided buffer,
        /// returns the number of bytes read, otherwise, returns 'null'.
        pub fn readNoEof(reader: anytype, buf: []u8) @TypeOf(reader).Error!?usize {
            const amt_read = try reader.readAll(buf);
            return if (amt_read < buf.len) amt_read else null;
        }

        /// Tries to read a BE integer from the stream; if the stream ends before supplying enough bytes
        /// for such an integer, returns null.
        pub fn readIntBigOrNull(reader: anytype, comptime T: type) @TypeOf(reader).Error!?T {
            var buffer: [(@typeInfo(T).Int.bits + 7) / 8]u8 = undefined;
            return if ((try readNoEof(reader, &buffer)) != null)
                null
            else
                std.mem.readIntBig(T, &buffer);
        }
    };
};
