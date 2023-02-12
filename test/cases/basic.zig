const std = @import("std");
const png = @import("../cases.zig").png;

const basic_png = @embedFile("basic.png");

test {
    const stderr = std.io.getStdErr().writer();
    try stderr.writeAll("\n");

    var fbs = std.io.fixedBufferStream(basic_png);
    var iter = png.chunkIterator(fbs.reader(), true);
    const data_reader = iter.dataReader();

    {
        var signature_check: png.ChunkIteratorSignatureCheck = undefined;
        iter.checkSignatureAdvanced(&signature_check) catch |err| switch (err) {}; // handle read errors
        try signature_check.unwrap();
    }

    while (true) {
        const header: png.ChunkHeader = blk: {
            var next_result: png.ChunkIteratorNextResult = undefined;
            iter.nextAdvanced(&next_result) catch |err| switch (err) {}; // handle read errors
            const maybe_header = try next_result.unwrap();
            break :blk maybe_header orelse break;
        };

        const data_buf = try std.testing.allocator.alloc(u8, header.length);
        defer std.testing.allocator.free(data_buf);

        try std.testing.expectEqual(data_buf.len, try data_reader.readAll(data_buf));
        try std.testing.expectEqual(@as(usize, 0), (try data_reader.readBoundedBytes(1)).len);

        // after all the chunk data has been read from
        // the data_reader, fetch the CRC code.
        if (!try iter.fetchExpectedCrc()) return error.FailedToFetchExpectedCrc;

        const expected_crc = iter.getExpectedCrc();
        const actual_crc = iter.getActualCrc();
        try std.testing.expectEqual(expected_crc, actual_crc);

        try stderr.print(
            \\{s}, length: {d}
            \\{any}
            \\CRC (expected): {any}
            \\CRC (actual):   {any}
            \\
            \\
        , .{
            header.type.string(),
            header.length,
            data_buf,
            expected_crc,
            actual_crc,
        });

        if (header.type == .IEND) break;
    }
}
