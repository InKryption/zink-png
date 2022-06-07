const std = @import("std");
const builtin = @import("builtin");

pub const png_signature: [8]u8 = .{ 137, 80, 78, 71, 13, 10, 26, 10 };

const ReadNoEofResult = enum { eof, no_eof };
fn readNoEof(reader: anytype, buf: []u8) @TypeOf(reader).Error!ReadNoEofResult {
    const amt_read = try reader.readAll(buf);
    return if (amt_read < buf.len) .eof else .no_eof;
}

fn skipBytesBuffered(reader: anytype, num_bytes: u64, buffer: []u8) @TypeOf(reader).Error!ReadNoEofResult {
    std.debug.assert(buffer.len >= 1);
    var remaining = num_bytes;

    while (remaining > 0) {
        const amt = std.math.min(remaining, buffer.len);
        switch (try readNoEof(reader, buffer[0..amt])) {
            .eof => return .eof,
            .no_eof => {},
        }
        remaining -= amt;
    }
    return .no_eof;
}

/// Tries to read a BE integer from the stream; if the stream ends before supplying enough bytes
/// for such an integer, returns null.
fn readIntBigOrNull(reader: anytype, comptime T: type) @TypeOf(reader).Error!?T {
    var buffer: [(@typeInfo(T).Int.bits + 7) / 8]u8 = undefined;
    return switch (try readNoEof(reader, &buffer)) {
        .eof => null,
        else => std.mem.readIntBig(T, &buffer),
    };
}

pub fn pngType(
    bytes: *const [4]u8,
) u32 {
    return std.mem.readIntBig(u32, bytes);
}
pub fn pngTypeStr(bytes: u32) [4]u8 {
    var result: [4]u8 = undefined;
    std.mem.writeIntBig(u32, &result, bytes);
    return result;
}

pub const PngRawChunkHeader = struct {
    length: u32,
    type: u32,
};
pub const PngRawChunk = struct {
    header: PngRawChunkHeader,
    p_data: [*]const u8,
    crc: u32,

    pub fn data(self: PngRawChunk) []const u8 {
        return self.p_data[0..self.header.length];
    }

    pub fn deinit(self: PngRawChunk, allocator: std.mem.Allocator) void {
        allocator.free(self.data());
    }

    const FmtConfig = struct {
        indentation: []const u8 = "    ",
        newline: []const u8 = "\n",
        data_width: u16 = 32,
    };

    pub fn fmt(
        self: PngRawChunk,
        config: FmtConfig,
    ) PngRawChunk.Formatter {
        return .{
            .value = self,
            .config = config,
        };
    }

    const Formatter = struct {
        value: PngRawChunk,
        config: FmtConfig,

        pub fn format(
            self: Formatter,
            comptime fmt_str: []const u8,
            options: std.fmt.FormatOptions,
            writer: anytype,
        ) @TypeOf(writer).Error!void {
            _ = fmt_str;
            _ = options;
            try writer.print(
                \\{{{0s}{1s}.type = '{2s}',{0s}{1s}.crc = 0x{3X},{0s}{1s}.data = [{4d}]u8{{
            , .{
                self.config.newline,
                self.config.indentation,
                pngTypeStr(self.value.header.type),
                self.value.crc,
                self.value.header.length,
            });
            for (self.value.data()) |byte, i| {
                if (i % self.config.data_width == 0) {
                    try writer.print("{0s}{1s}{1s}", .{
                        self.config.newline,
                        self.config.indentation,
                    });
                }
                try writer.print("{d}, ", .{byte});
            }
            if (self.value.header.length != 0) {
                try writer.print("{s}{s}", .{ self.config.newline, self.config.indentation });
            }
            try writer.print("}},{s}}}", .{self.config.newline});
        }
    };
};

pub const PngRawChunkStreamState = enum {
    begin,
    reading,
    end,
};
pub const DataCapture = union(enum) {
    /// if used, the caller should `deinit` the returned chunk.
    allocator: std.mem.Allocator,
    /// will be used to skip over the data.
    dont_capture: []u8,
};
pub fn PngRawChunkStream(comptime ReaderType: type) type {
    return struct {
        const Self = @This();
        reader: ReaderType,
        state: PngRawChunkStreamState,

        pub fn init(reader: ReaderType) Self {
            return Self{
                .reader = reader,
                .state = .begin,
            };
        }

        pub const Error = ReaderType.Error || std.mem.Allocator.Error || error{
            StreamTooShortForPngSignature,
            IncorrectPngSignature,
            StreamEndedMidChunk,
            SuppliedBufferTooShortForChunkData,
        };
        pub fn nextAlloc(self: *Self, allocator: std.mem.Allocator) ?(Error!PngRawChunk) {
            return self.nextAdvanced(.{ .allocator = allocator });
        }
        pub fn nextSkip(self: *Self, options: ReaderType.SkipBytesOptions) ?(Error!PngRawChunk) {
            var buf: [options.buf_size]u8 = undefined;
            return self.nextAdvanced(.{ .dont_capture = &buf });
        }

        pub fn nextAdvanced(
            self: *Self,
            data_capture: DataCapture,
        ) ?(Error!PngRawChunk) {
            // mark stream as done if an error is encountered.
            errdefer self.state = .end;

            switch (self.state) {
                .begin => {
                    var actual: [png_signature.len]u8 = undefined;

                    switch (try readNoEof(self.reader, &actual)) {
                        .eof => return error.StreamTooShortForPngSignature,
                        .no_eof => {},
                    }

                    if (!std.mem.eql(u8, &actual, &png_signature)) {
                        return error.IncorrectPngSignature;
                    }

                    self.state = .reading;
                },
                .reading => {},
                .end => return null,
            }

            const header = PngRawChunkHeader{
                .length = (try readIntBigOrNull(self.reader, u32)) orelse {
                    self.state = .end;
                    return null;
                },
                .type = (try readIntBigOrNull(self.reader, u32)) orelse
                    return error.StreamEndedMidChunk,
            };
            const data: [*]const u8 = switch (data_capture) {
                .allocator => |allocator| blk: {
                    const data = try allocator.alloc(u8, header.length);
                    errdefer allocator.free(data);

                    switch (try readNoEof(self.reader, data)) {
                        .eof => return error.StreamEndedMidChunk,
                        .no_eof => {},
                    }

                    break :blk data.ptr;
                },
                .dont_capture => |buf| blk: {
                    switch (try skipBytesBuffered(self.reader, header.length, buf)) {
                        .eof => return error.StreamEndedMidChunk,
                        .no_eof => {},
                    }
                    break :blk undefined;
                },
            };
            const crc = (try readIntBigOrNull(self.reader, u32)) orelse return error.StreamEndedMidChunk;

            return PngRawChunk{
                .header = header,
                .p_data = data,
                .crc = crc,
            };
        }
    };
}
pub fn pngRawChunkStream(reader: anytype) PngRawChunkStream(@TypeOf(reader)) {
    return PngRawChunkStream(@TypeOf(reader)).init(
        reader,
    );
}

test {
    std.debug.print("\n", .{});
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const data = @embedFile("img1.png");
    var data_stream = std.io.fixedBufferStream(data);

    var chunk_stream = pngRawChunkStream(data_stream.reader());
    while (chunk_stream.nextAlloc(arena.allocator())) |maybe_chunk| {
        const chunk = try maybe_chunk;
        std.debug.print("{any}\n", .{chunk.fmt(.{})});
    }
}
