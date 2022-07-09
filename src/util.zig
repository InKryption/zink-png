const std = @import("std");

pub fn ErrorSetFromValue(comptime value: anyerror) type {
    const set = [_]std.builtin.Type.Error{std.builtin.Type.Error{ .name = @errorName(value) }};
    return @Type(.{ .ErrorSet = &set });
}

pub fn NormalizedErrorSet(comptime ErrorSet: type) type {
    const info: []const std.builtin.Type.Error = @typeInfo(ErrorSet).ErrorSet orelse return anyerror;
    var values: [info.len]anyerror = .{undefined} ** info.len;
    for (values) |*val, i| val.* = @field(anyerror, info[i].name);
    std.sort.sort(anyerror, &values, void{}, struct {
        fn lessThan(ctx: void, lhs: anyerror, rhs: anyerror) bool {
            ctx;
            return std.mem.lessThan(u8, @errorName(lhs), @errorName(rhs));
        }
    }.lessThan);
    return NormalizedErrorSetImpl(values.len, values);
}

inline fn NormalizedErrorSetImpl(
    comptime error_value_count: comptime_int,
    comptime error_values: [error_value_count]anyerror,
) type {
    var ErrorSet = error{};
    for (error_values) |val| {
        ErrorSet = ErrorSet || ErrorSetFromValue(val);
    }
    return ErrorSet;
}

pub fn as(comptime T: type, value: anytype) T {
    return value;
}

pub const io = struct {
    pub fn errorLimitedReader(inner: anytype, limit: u64) ErrorLimitedReader(@TypeOf(inner)) {
        return ErrorLimitedReader(@TypeOf(inner)).init(inner, limit);
    }
    pub fn ErrorLimitedReader(comptime InnerReader: type) type {
        return struct {
            const Self = @This();
            inner: std.io.CountingReader(InnerReader),
            limit: u64,

            pub fn init(inner: InnerReader, limit: u64) Self {
                return Self{
                    .inner = std.io.countingReader(inner),
                    .limit = limit,
                };
            }

            pub const ReadError = error{AttemptedToReadMoreThanLimit} || InnerReader.Error;
            pub const Reader = std.io.Reader(*Self, ReadError, Self.read);
            pub fn reader(self: *Self) Self.Reader {
                return .{ .context = self };
            }

            fn read(self: *Self, buffer: []u8) ReadError!usize {
                std.debug.assert(self.inner.bytes_read <= self.limit);
                if (self.inner.bytes_read == self.limit) {
                    return error.AttemptedToReadMoreThanLimit;
                }
                const remaining_bytes = self.limit - self.inner.bytes_read;
                return self.inner.reader().read(buffer[0..std.math.min(remaining_bytes, buffer.len)]);
            }
        };
    }

    pub fn readIntoBoundedArray(
        reader: anytype,
        /// Must be `std.BoundedArray(u8, n)`.
        bounded_array: anytype,
    ) @TypeOf(reader).Error!void {
        if (*std.BoundedArray(u8, bounded_array.buffer.len) != @TypeOf(bounded_array)) {
            @compileError("Expected `*std.BoundedArray(u8, n)`, got `" ++ @typeName(bounded_array) ++ "`.");
        }

        while (bounded_array.len != bounded_array.buffer.len) {
            const amt = try reader.read(bounded_array.unusedCapacitySlice());
            if (amt == 0) break;
            bounded_array.resize(bounded_array.len + amt) catch unreachable;
        }
    }
    pub fn readBoundedArray(reader: anytype, comptime max_bytes: usize) @TypeOf(reader).Error!std.BoundedArray(u8, max_bytes) {
        var result = std.BoundedArray(u8, max_bytes){};
        try readIntoBoundedArray(reader, &result);
        return result;
    }

    pub fn readBytesNoEof(reader: anytype, comptime num_bytes: usize) @TypeOf(reader).Error!?[num_bytes]u8 {
        var buffer: [num_bytes]u8 = undefined;
        const bytes_read = try reader.readAll(buffer[0..]);
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
