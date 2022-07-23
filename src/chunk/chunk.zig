pub const Header = @import("Header.zig");
pub const IHDR = @import("IHDR.zig");
pub const PLTE = @import("PLTE.zig");

comptime {
    _ = Header;
    _ = IHDR;
    _ = PLTE;
}
