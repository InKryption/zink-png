const std = @import("std");
const png = @import("../cases.zig").png;

const basic_png = @embedFile("basic.png");

test {
    var fbs = std.io.fixedBufferStream(basic_png);
    var iter = png.chunkIterator(fbs.reader());
    const data_reader = iter.reader();

    try (try iter.checkSignature()).unwrap();
    while (try (try iter.next()).unwrap()) |header| {
        try data_reader.skipBytes(header.length, .{});
        if (!try iter.fetchExpectedCrc()) return error.FailedToFetchExpectedCrc;
        const crc = iter.getCrc();
        try std.testing.expectEqual(crc.expected, crc.actual);
        if (header.type == .IEND) break;
    }
}
