const std = @import("std");
const png = @import("../cases.zig").png;

const basic_png = @embedFile("basic.png");

fn expectChunkData(
    chunk_iter: anytype,
    expected_header: png.ChunkHeader,
    expected_data: []const u8,
) !void {
    const data_reader = chunk_iter.dataReader();

    const actual_header: png.ChunkHeader = blk: {
        var next_result: png.ChunkIteratorNextResult = undefined;
        chunk_iter.nextAdvanced(&next_result) catch |err| switch (err) {}; // handle read errors
        const maybe_header = try next_result.unwrap();
        break :blk maybe_header orelse return;
    };
    try std.testing.expectEqual(expected_header, actual_header);
    try std.testing.expect(try data_reader.isBytes(expected_data));
    try std.testing.expectEqual(@as(usize, 0), (try data_reader.readBoundedBytes(1)).len);

    if (!try chunk_iter.fetchExpectedCrc()) return error.FailedToFetchExpectedCrc;
    const expected_crc = chunk_iter.getExpectedCrc();
    const actual_crc = chunk_iter.getActualCrc();
    try std.testing.expectEqual(expected_crc, actual_crc);
}

test {
    var fbs = std.io.fixedBufferStream(basic_png);
    var iter = png.chunkIterator(fbs.reader(), true);

    try (try iter.checkSignature()).unwrap();
    while (true) {
        const header: png.ChunkHeader = blk: {
            var next_result: png.ChunkIteratorNextResult = undefined;
            iter.nextAdvanced(&next_result) catch |err| switch (err) {}; // handle read errors
            const maybe_header = try next_result.unwrap();
            break :blk maybe_header orelse break;
        };

        const chunk_data = try std.testing.allocator.alloc(u8, header.length);
        defer std.testing.allocator.free(chunk_data);

        try std.testing.expectEqual(chunk_data.len, try iter.dataReader().readAll(chunk_data));
        try std.testing.expectEqual(@as(usize, 0), (try iter.dataReader().readBoundedBytes(1)).len);

        // after all the chunk data has been read from
        // the data_reader, fetch the CRC code.
        if (!try iter.fetchExpectedCrc()) return error.FailedToFetchExpectedCrc;

        const real_expected_crc = std.hash.Crc32.hash(chunk_data);
        const expected_crc = iter.getExpectedCrc();
        const actual_crc = iter.getActualCrc();

        try std.testing.expectEqual(real_expected_crc, expected_crc);
        try std.testing.expectEqual(expected_crc, actual_crc);

        if (header.type == .IEND) break;
    }

    // test contents

    fbs.reset();
    iter = png.chunkIterator(fbs.reader(), true);

    {
        var signature_check: png.ChunkIteratorSignatureCheck = undefined;
        iter.checkSignatureAdvanced(&signature_check) catch |err| switch (err) {}; // handle read errors
        try signature_check.unwrap();
    }

    try expectChunkData(
        &iter,
        .{ .type = .IHDR, .length = 13 },
        &png.IHDR.Data.toBytes(.{
            .width = 10,
            .height = 10,
            .bit_depth = .@"8",
            .color_type = .grayscale,
            .compression_method = .method_0,
            .filter_method = .method_0,
            .interlace_method = .none,
        }),
    );
    try expectChunkData(
        &iter,
        .{ .type = .IDAT, .length = 121 },
        &[121]u8{
            8,   215, 1,   110, 0,   145, 255, 0,   243, 0,   177,
            56,  145, 81,  210, 0,   146, 42,  2,   174, 20,  14,
            24,  91,  23,  205, 24,  15,  234, 2,   7,   241, 32,
            4,   173, 160, 5,   25,  14,  52,  3,   125, 231, 219,
            2,   233, 46,  48,  235, 240, 226, 3,   20,  241, 236,
            241, 137, 202, 206, 243, 248, 226, 2,   221, 226, 113,
            11,  149, 254, 253, 245, 58,  1,   3,   171, 203, 224,
            17,  0,   1,   45,  228, 6,   218, 3,   16,  255, 254,
            252, 103, 173, 236, 235, 255, 196, 4,   237, 228, 144,
            70,  24,  3,   4,   245, 144, 214, 4,   133, 19,  43,
            44,  227, 3,   158, 1,   165, 239, 125, 121, 51,  177,
        },
    );
    try expectChunkData(
        &iter,
        .{ .type = .IEND, .length = 0 },
        &[0]u8{},
    );
}
