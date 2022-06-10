const std = @import("std");

pub fn ErrorSetFromValue(comptime value: anyerror) type {
    const set = [_]std.builtin.Type.Error{std.builtin.Type.Error{ .name = @errorName(value) }};
    return @Type(.{ .ErrorSet = &set });
}

pub fn MemoizedErrorSet(comptime ErrorSet: type) type {
    const info: []const std.builtin.Type.Error = @typeInfo(ErrorSet).ErrorSet orelse return anyerror;
    var values: [info.len]anyerror = .{undefined} ** info.len;
    for (values) |*val, i| val.* = @field(anyerror, info[i].name);
    std.sort.sort(anyerror, &values, void{}, struct {
        fn lessThan(ctx: void, lhs: anyerror, rhs: anyerror) bool {
            ctx;
            return std.mem.lessThan(u8, @errorName(lhs), @errorName(rhs));
        }
    }.lessThan);
    return MemoizedErrorSetImpl(values.len, values);
}

inline fn MemoizedErrorSetImpl(
    comptime error_value_count: comptime_int,
    comptime error_values: [error_value_count]anyerror,
) type {
    var ErrorSet = error{};
    for (error_values) |val| {
        ErrorSet = ErrorSet || ErrorSetFromValue(val);
    }
    return ErrorSet;
}

pub const io = struct {
    pub fn ReadAllResult(comptime ErrorSet: type) type {
        return ReadAllResultImpl(MemoizedErrorSet(ErrorSet));
    }
    fn ReadAllResultImpl(comptime ErrorSet: type) type {
        return struct {
            bytes_read: usize,
            err: if (@sizeOf(Err) == 0) Err else ?Err,

            pub inline fn getError(self: @This()) ?Err {
                comptime if (@sizeOf(Err) == 0) return null;
                return self.err;
            }

            pub const Err = ErrorSet;
            pub fn unwrap(self: @This()) Err!usize {
                return self.getError() orelse self.bytes_read;
            }
        };
    }
    /// Returns the number of bytes read + the error, if any. If the number read is smaller than `buffer.len`, it
    /// means the stream reached the end. Reaching the end of a stream is not an error
    /// condition.
    /// If an error did occur, the number of bytes read will be less than `buffer.len`.
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
        const bytes_read = try readAll(reader, &buffer).unwrap();
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
