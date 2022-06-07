const std = @import("std");
const builtin = @import("builtin");

pub fn validPngSignature(reader: anytype) bool {
    const png_signature: [8]u8 = .{ 137, 80, 78, 71, 13, 10, 26, 10 };
    return reader.isBytes(&png_signature) catch false;
}
pub fn pngType(
    byte_a: u8,
    byte_b: u8,
    byte_c: u8,
    byte_d: u8,
) u32 {
    return std.mem.readIntBig(u32, &[4]u8{ byte_a, byte_b, byte_c, byte_d });
}
pub fn pngTypeStr(bytes: u32) [4]u8 {
    var result: [4]u8 = undefined;
    std.mem.writeIntBig(u32, &result, bytes);
    return result;
}

pub const PngChunkHeader = struct {
    length: u32,
    type: u32,
};

pub const ReadPngChunkHeaderError = error{
    PngChunkHeaderLengthBytesUnavailable,
    PngChunkHeaderTypeBytesUnavailable,
};
pub fn readPngChunkHeader(reader: anytype) (@TypeOf(reader).Error || ReadPngChunkHeaderError)!PngChunkHeader {
    return PngChunkHeader{
        .length = reader.readIntBig(u32) catch return error.PngChunkHeaderLengthBytesUnavailable,
        .type = reader.readIntBig(u32) catch return error.PngChunkHeaderTypeBytesUnavailable,
    };
}

pub const ReadPngChunkCrcError = error{
    PngChunkCrcUnavailable,
};
pub fn readPngChunkCrc(reader: anytype) (@TypeOf(reader).Error || ReadPngChunkCrcError)!u32 {
    return reader.readIntBig(u32) catch |err| switch (err) {
        error.EndOfStream => error.PngChunkCrcUnavailable,
        else => |e| return e,
    };
}

pub const PngRawChunk = struct {
    header: PngChunkHeader,
    p_data: [*]const u8,
    crc: u32,

    pub fn deinit(self: PngRawChunk, allocator: std.mem.Allocator) void {
        allocator.free(self.data());
    }

    pub fn data(self: PngRawChunk) []const u8 {
        return self.p_data[0..self.header.length];
    }
};

pub const PngRawChunkStreamState = enum {
    begin,
    reading,
    end,
};
pub const DataCapture = union(enum) {
    /// if used, the caller should `deinit` the returned chunk.
    allocator: std.mem.Allocator,
    buffer: []u8,
    assert_empty,
    dont_capture,
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

        pub const Error = ReaderType.Error || ReadPngChunkHeaderError || ReadPngChunkCrcError || std.mem.Allocator.Error || error{
            BadPngSignature,
            StreamEndedMidChunk,
            SuppliedBufferTooShortForChunkData,
            AssertedNonEmptyDataAsEmpty,
        };
        pub fn next(
            self: *Self,
            data_capture: DataCapture,
            comptime skip_bytes_config: ReaderType.SkipBytesOptions,
        ) ?(Error!PngRawChunk) {
            switch (self.state) {
                .begin => {
                    if (!validPngSignature(self.reader)) {
                        self.state = .end;
                        return error.BadPngSignature;
                    }
                    self.state = .reading;
                },
                .reading => {},
                .end => return null,
            }

            const header = try readPngChunkHeader(self.reader);
            const data: [*]const u8 = switch (data_capture) {
                .allocator => |allocator| blk: {
                    const data = try allocator.alloc(u8, header.length);
                    errdefer allocator.free(data);

                    self.reader.readNoEof(data) catch |err| return switch (err) {
                        error.EndOfStream => return error.StreamEndedMidChunk,
                        else => |e| e,
                    };

                    break :blk data.ptr;
                },
                .buffer => |buffer| blk: {
                    if (buffer.len < header.length) {
                        return error.SuppliedBufferTooShortForChunkData;
                    }
                    self.reader.readNoEof(buffer[0..header.length]) catch |err| return switch (err) {
                        error.EndOfStream => return error.StreamEndedMidChunk,
                        else => |e| e,
                    };
                    break :blk buffer.ptr;
                },
                .assert_empty => blk: {
                    if (header.length != 0) {
                        return error.AssertedNonEmptyDataAsEmpty;
                    }
                    break :blk undefined;
                },
                .dont_capture => blk: {
                    self.reader.skipBytes(header.length, skip_bytes_config) catch |err| return switch (err) {
                        error.EndOfStream => return error.StreamEndedMidChunk,
                        else => |e| e,
                    };
                    break :blk undefined;
                },
            };
            const crc = try readPngChunkCrc(self.reader);

            switch (header.type) {
                pngType('I', 'E', 'N', 'D') => {
                    self.state = .end;
                },
                else => {},
            }

            return PngRawChunk{
                .header = header,
                .p_data = data,
                .crc = crc,
            };
        }
    };
}
pub fn pngRawChunkStream(reader: anytype) PngRawChunkStream(@TypeOf(reader)) {
    return PngRawChunkStream(@TypeOf(reader)).init(reader);
}

test {}
