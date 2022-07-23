const std = @import("std");
const png = @import("../main.zig");
const chunk = png.chunk;

// -- Start ColorType.palette (3) --

/// Represents the alpha channels corresponding
/// to the palette entries stored in the tRNS chunk when an image
/// has a color type 3 (png.chunk.IHDR.ColorType.palette).
pub const AlphaChannelEntries = std.BoundedArray(u8, 256);

/// Returns null if `plte.count() < plte_alphas.len`. Otherwise returns
/// a bounded array consisting of the alpha channel values present in `plte_alphas`, followed
/// by `plte.count() - plte_alphas.len` alpha channel values, all equal to `std.math.maxInt(u8)`.
pub fn fullCorrespondingAlphaChannels(plte: chunk.PLTE.ChunkData, plte_alphas: AlphaChannelEntries) ?AlphaChannelEntries {
    if (plte.count() < plte_alphas.len) return null;

    var result: AlphaChannelEntries = plte_alphas;
    result.appendNTimesAssumeCapacity(std.math.maxInt(u8), plte.count() - plte_alphas.len);
    return result;
}

test "tRNS - color type 3" {
    _ = AlphaChannelEntries;
    _ = fullCorrespondingAlphaChannels;
    std.log.warn("TODO: Add actual testing here", .{});
    return error.SkipZigTest;
}

// -- End ColorType.palette (3) --

/// Returns an integer that can hold values starting from 0, and up to `(2^bit_depth)-1`.
pub fn BitDepthRangedInt(comptime bit_depth: chunk.IHDR.BitDepth) type {
    return std.math.IntFittingRange(0, std.math.ceilPowerOfTwoPromote(u8, @enumToInt(bit_depth)) - 1);
}

// -- Start ColorType.grayscale --

/// Represents the gray level value stored in the tRNS chunk
/// when an image has a color type 0 (png.chunk.IHDR.ColorType.grayscale).
/// The tag corresponds to the bit depth given in the IHDR chunk.
pub const GrayLevel = union(chunk.IHDR.BitDepth) {
    bd1: BitDepthRangedInt(.bd1),
    bd2: BitDepthRangedInt(.bd2),
    bd4: BitDepthRangedInt(.bd4),
    bd8: BitDepthRangedInt(.bd8),
    bd16: BitDepthRangedInt(.bd16),

    /// Returns null if the integer would not fit into an integer `BitDepthRangedInt(bit_depth)`.
    pub fn fromBytesWithBitDepth(bytes: [2]u8, bit_depth: chunk.IHDR.BitDepth) ?GrayLevel {
        const raw_value = std.mem.readIntBig(u16, bytes[0..]);
        return switch (bit_depth) {
            .bd1 => GrayLevel{ .bd1 = std.math.cast(BitDepthRangedInt(.bd1), raw_value) orelse return null },
            .bd2 => GrayLevel{ .bd2 = std.math.cast(BitDepthRangedInt(.bd2), raw_value) orelse return null },
            .bd4 => GrayLevel{ .bd4 = std.math.cast(BitDepthRangedInt(.bd4), raw_value) orelse return null },
            .bd8 => GrayLevel{ .bd8 = std.math.cast(BitDepthRangedInt(.bd8), raw_value) orelse return null },
            .bd16 => GrayLevel{ .bd16 = @intCast(BitDepthRangedInt(.bd16), raw_value) },
        };
    }

    pub fn toBytes(gray_level: GrayLevel) [2]u8 {
        var result: [2]u8 align(@alignOf(u16)) = .{0} ** 2;
        switch (gray_level) {
            .bd1 => |value| std.mem.writeIntBig(u16, &result, value),
            .bd2 => |value| std.mem.writeIntBig(u16, &result, value),
            .bd4 => |value| std.mem.writeIntBig(u16, &result, value),
            .bd8 => |value| std.mem.writeIntBig(u16, &result, value),
            .bd16 => |value| std.mem.writeIntBig(u16, &result, value),
        }
        return result;
    }
};

test "tRNS - color type 0" {
    _ = GrayLevel;
    std.log.warn("TODO: Add actual testing here", .{});
    return error.SkipZigTest;
}

// -- End ColorType.grayscale --

// -- Start ColorType.truecolor (2) --

/// Represents the rgb level values stored in the tRNS chunk
/// when an image has a color type 2 (png.chunk.IHDR.ColorType.truecolor).
/// The tag corresponds to the bit depth given in the IHDR chunk.
pub const ColorLevels = union(chunk.IHDR.BitDepth) {
    bd1: Rgb(.bd1),
    bd2: Rgb(.bd2),
    bd4: Rgb(.bd4),
    bd8: Rgb(.bd8),
    bd16: Rgb(.bd16),

    /// Returns an array of integers, wherein:
    /// r: array[0]
    /// g: array[1]
    /// b: array[2]
    pub fn Rgb(comptime bit_depth: chunk.IHDR.BitDepth) type {
        return [3]BitDepthRangedInt(bit_depth);
    }

    pub const FromBytesWithDepthError = error{ RedOverflow, GreenOverflow, BlueOverflow };
    pub fn fromBytesWithBitDepth(bytes: [3 * 2]u8, bit_depth: chunk.IHDR.BitDepth) FromBytesWithDepthError!ColorLevels {
        const raw_values: [3]u16 = raw_values: {
            var raw_values = std.mem.bytesToValue([3]u16, bytes);
            for (raw_values) |*rv| rv.* = std.mem.bigToNative(u16, rv.*);
            break :raw_values raw_values;
        };
        return switch (bit_depth) {
            .bd1 => ColorLevels{ .bd1 = try castArray(.bd1, raw_values) },
            .bd2 => ColorLevels{ .bd2 = try castArray(.bd2, raw_values) },
            .bd4 => ColorLevels{ .bd4 = try castArray(.bd4, raw_values) },
            .bd8 => ColorLevels{ .bd8 = try castArray(.bd8, raw_values) },
            .bd16 => ColorLevels{ .bd16 = raw_values },
        };
    }

    pub fn toBytes(color_levels: ColorLevels) [3 * 2]u8 {
        return switch (color_levels) {
            .bd1 => |rgb| makeBigIntArray(.bd1, rgb),
            .bd2 => |rgb| makeBigIntArray(.bd2, rgb),
            .bd4 => |rgb| makeBigIntArray(.bd4, rgb),
            .bd8 => |rgb| makeBigIntArray(.bd8, rgb),
            .bd16 => |rgb| makeBigIntArray(.bd16, rgb),
        };
    }

    fn castArray(
        comptime bit_depth: chunk.IHDR.BitDepth,
        raw_values: [3]u16,
    ) FromBytesWithDepthError![3]BitDepthRangedInt(bit_depth) {
        const DstInt = BitDepthRangedInt(bit_depth);
        return [3]DstInt{
            std.math.cast(DstInt, raw_values[0]) orelse return error.RedOverflow,
            std.math.cast(DstInt, raw_values[1]) orelse return error.GreenOverflow,
            std.math.cast(DstInt, raw_values[2]) orelse return error.BlueOverflow,
        };
    }
    fn makeBigIntArray(
        comptime bit_depth: chunk.IHDR.BitDepth,
        values: Rgb(bit_depth),
    ) [3 * 2]u8 {
        var result: [3 * 2]u8 = .{0} ** 6;

        var result_fbs = std.io.fixedBufferStream(result[0..]);
        defer std.debug.assert(result_fbs.pos == result.len);

        inline for (values) |val| {
            result_fbs.writer().writeIntBig(u16, val) catch unreachable;
        }

        return result;
    }
};

test "tRNS - color type 2" {
    _ = ColorLevels;
    std.log.warn("TODO: Add actual testing here", .{});
    return error.SkipZigTest;
}

// -- End ColorType.truecolor (2) --
