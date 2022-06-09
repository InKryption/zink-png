const std = @import("std");
const png = @import("main.zig");

const ChunkHeader = png.ChunkHeader;
const ChunkType = png.ChunkType;

pub const RawChunk = struct {
    header: ChunkHeader,
    p_data: [*]const u8,
    crc: u32,

    pub fn data(self: RawChunk) []const u8 {
        return self.p_data[0..self.header.length];
    }

    /// Should only be called if returned from a call to `RawChunkStream.next`, where an allocator was supplied;
    pub fn deinit(self: RawChunk, allocator: std.mem.Allocator) void {
        allocator.free(self.data());
    }
};

pub fn rawChunkStream(reader: anytype) RawChunkStream(@TypeOf(reader)) {
    return RawChunkStream(@TypeOf(reader)).init(reader);
}

pub fn RawChunkStream(comptime ReaderType: type) type {
    return struct {
        const Self = @This();
        reader: Reader,
        state: RawChunkStreamState,

        pub const Reader = ReaderType;
        pub const State = RawChunkStreamState;
        pub const StartResult = RawChunkStreamStartResult(ReaderType.Error);
        pub const NextResult = RawChunkStreamNextResult(ReaderType.Error);

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

                    var bytes: [png.signature.len]u8 = undefined;
                    const bytes_read = self.reader.readAll(&bytes) catch |err| {
                        return StartResult{ .initial_read_fail = .{ .err = err } };
                    };
                    if (bytes_read != bytes.len) {
                        return StartResult{ .initial_read_fail = .{ .err = null } };
                    }

                    if (!std.mem.eql(u8, &bytes, &png.signature)) {
                        return .bad_png_signature;
                    }

                    self.state = .in_progress;
                },
                .in_progress => unreachable,
                .end => unreachable,
            }
            return .ok;
        }

        pub const MaybeAllocator = union(enum) {
            allocator: std.mem.Allocator,
            skip: []u8,
        };

        /// Returns 'null' when the stream ends after stepping through only whole chunks,
        /// or after failing to read a whole chunk.
        /// If no allocator is supplied, on success, the returned chunk will have an undefined pointer,
        /// and thus should not be used.
        /// Must have called `start` beforehand.
        pub fn next(self: *Self, maybe_allocator: MaybeAllocator) ?NextResult {
            switch (self.state) {
                .start => unreachable,
                .in_progress => {},
                .end => return null,
            }

            // Set state to 'end' by default; if nothing impedes control flow from reaching the end of the function,
            // it will be set back to `in_progress`.
            self.state = .end;

            const header: ChunkHeader = header: {
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

                break :header ChunkHeader{
                    .length = length,
                    .type = @intToEnum(ChunkType, @"type"),
                };
            };

            const data: []const u8 = switch (maybe_allocator) {
                .allocator => |allocator| data: {
                    const data: []u8 = allocator.alloc(u8, header.length) catch |err| switch (err) {
                        error.OutOfMemory => return NextResult{ .out_of_mem_for_data = NextResult.OutOfMemForData{ .header = header } },
                    };
                    const read_all_result = util.io.readAll(self.reader, data);
                    if (read_all_result.bytes_read < data.len) {
                        return NextResult{ .partial_data_bytes = NextResult.PartialDataBytes{
                            .header = header,
                            .bytes = allocator.shrink(data, read_all_result.bytes_read),
                            .err = read_all_result.err,
                        } };
                    }
                    std.debug.assert(read_all_result.err == null);

                    break :data data;
                },
                .skip => |buffer| data: {
                    var remaining: usize = header.length;
                    while (remaining > 0) {
                        const amt = std.math.min(remaining, buffer.len);
                        const read_all_result = util.io.readAll(self.reader, buffer[0..amt]);

                        if (read_all_result.bytes_read < amt) {
                            remaining -= read_all_result.bytes_read;
                            return NextResult{
                                .partial_data_bytes_no_capture = NextResult.PartialDataBytesNoCapture{
                                    .header = header,
                                    .bytes_len = header.length - remaining,
                                    .err = read_all_result.err,
                                },
                            };
                        }
                        std.debug.assert(read_all_result.err == null);
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
                .ok = RawChunk{
                    .header = header,
                    .p_data = data.ptr,
                    .crc = crc,
                },
            };
        }
    };
}

const RawChunkStreamState = enum {
    start,
    in_progress,
    end,
};
const RawChunkStreamStartResultTag = enum {
    ok,
    initial_read_fail,
    bad_png_signature,
};
const RawChunkStreamNextResultTag = enum {
    ok,
    no_length_bytes,
    no_type_bytes,
    out_of_mem_for_data,
    partial_data_bytes_no_capture,
    partial_data_bytes,
    no_crc_bytes_no_capture,
    no_crc_bytes,
};

fn RawChunkStreamStartResult(comptime ReadError: type) type {
    return union(RawChunkStreamStartResultTag) {
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
fn RawChunkStreamNextResult(comptime ReadError: type) type {
    return union(RawChunkStreamNextResultTag) {
        const Self = @This();
        ok: RawChunk,
        no_length_bytes: NoLengthBytes,
        no_type_bytes: NoTypeBytes,
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
        pub const OutOfMemForData = struct {
            header: ChunkHeader,
        };
        pub const PartialDataBytesNoCapture = struct {
            header: ChunkHeader,
            bytes_len: usize,
            /// 'null' if stream ended.
            err: ?Error,
        };
        pub const PartialDataBytes = struct {
            header: ChunkHeader,
            bytes: []const u8,
            /// 'null' if stream ended.
            err: ?Error,
        };
        pub const NoCrcBytes = struct {
            header: ChunkHeader,
            p_data: [*]const u8,
            /// 'null' if stream ended.
            err: ?Error,
        };
        pub const NoCrcBytesNoCapture = HeaderMaybeErrOnly;

        const HeaderMaybeErrOnly = struct {
            header: ChunkHeader,
            /// 'null' if stream ended.
            err: ?Error,
        };
    };
}

const util = struct {
    const io = struct {
        pub fn ReadAllResult(comptime ErrorSet: type) type {
            return struct {
                bytes_read: usize,
                err: ?ErrorSet,

                pub fn unwrap(self: @This()) ErrorSet!usize {
                    return self.err orelse self.bytes_read;
                }
            };
        }
        pub fn readAll(reader: anytype, buffer: []u8) ReadAllResult(@TypeOf(reader).Error) {
            const Result = ReadAllResult(@TypeOf(reader).Error);

            var index: usize = 0;
            while (index != buffer.len) {
                const amt = reader.read(buffer[index..]) catch |err|
                    return Result{ .bytes_read = index, .err = err };
                if (amt == 0) break;
                index += amt;
            }

            return Result{ .bytes_read = index, .err = null };
        }

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
