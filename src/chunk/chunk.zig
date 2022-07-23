const std = @import("std");

pub const Header = @import("Header.zig");
pub const IHDR = @import("IHDR.zig");
pub const PLTE = @import("PLTE.zig");

pub const Type = enum(u32) {
    pub const Tag = @typeInfo(Type).Enum.tag_type;
    // Standard Critical Chunk Types
    IHDR = std.mem.readIntBig(u32, "IHDR"),
    PLTE = std.mem.readIntBig(u32, "PLTE"),
    IDAT = std.mem.readIntBig(u32, "IDAT"),
    IEND = std.mem.readIntBig(u32, "IEND"),

    // Standard Ancillary Chunk Types
    bKGD = std.mem.readIntBig(u32, "bKGD"),
    cHRM = std.mem.readIntBig(u32, "cHRM"),
    gAMA = std.mem.readIntBig(u32, "gAMA"),
    hIST = std.mem.readIntBig(u32, "hIST"),
    pHYs = std.mem.readIntBig(u32, "pHYs"),
    sBIT = std.mem.readIntBig(u32, "sBIT"),
    tEXt = std.mem.readIntBig(u32, "tEXt"),
    tIME = std.mem.readIntBig(u32, "tIME"),
    tRNS = std.mem.readIntBig(u32, "tRNS"),
    zTXt = std.mem.readIntBig(u32, "zTXt"),

    // Non-Standard Chunk Types
    _,

    pub fn from(bytes: [4]u8) Type {
        return @intToEnum(Type, (std.mem.readIntBig(u32, bytes[0..])));
    }

    pub fn string(self: Type) [4]u8 {
        return std.mem.toBytes(self.intBig());
    }

    pub fn isValidAscii(self: Type) bool {
        for (self.string()) |byte| {
            if (!std.ascii.isAlpha(byte)) return false;
        }
        return true;
    }

    pub fn int(self: Type) Tag {
        return @enumToInt(self);
    }
    pub fn intBig(self: Type) Tag {
        return std.mem.nativeToBig(Tag, self.int());
    }
    pub fn intLittle(self: Type) Tag {
        return std.mem.nativeToLittle(Tag, self.int());
    }

    pub fn property(self: Type, byte_index: u2) bool {
        return (self.string()[byte_index] & 32) != 0;
    }

    pub fn format(
        ch_type: Type,
        comptime fmt_str: []const u8,
        fmt_options: std.fmt.FormatOptions,
        writer: anytype,
    ) @TypeOf(writer).Error!void {
        _ = fmt_options;
        const lazy = struct {
            const unknown_format_str = @compileError("Unknown format string '" ++ fmt_str ++ "'.\n");
        };

        switch (fmt_str.len) {
            0 => {},
            else => lazy.unknown_format_str,
        }
        return writer.print("{s}({s})", .{ @typeName(@This()), ch_type.string() });
    }
};

comptime {
    _ = Header;
    _ = Type;
    _ = IHDR;
    _ = PLTE;
}
