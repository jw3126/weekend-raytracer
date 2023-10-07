const std = @import("std");
const sqrt = std.math.sqrt;
const floatEps = std.math.floatEps;
const pow = std.math.pow;

pub fn SVec(comptime len: comptime_int, comptime T: type) type {
    return struct {
        const Self = @This();
        data: [len]T,

        fn undef() Self {
            var ret = Self{
                .data = undefined,
            };
            return ret;
        }

        pub fn fill(value: T) Self {
            var ret = Self.undef();
            @memset(&ret.data, value);
            return ret;
        }

        pub fn fromSlice(slice: []const T) Self {
            var ret = Self.undef();
            for (0..len) |i| {
                ret.data[i] = slice[i];
            }
            return ret;
        }

        pub fn at(self: Self, index: anytype) T {
            return self.data[index];
        }

        pub fn add(self: Self, other: Self) Self {
            var ret = Self.undef();
            for (0..len) |i| {
                ret.data[i] = self.data[i] + other.data[i];
            }
            return ret;
        }

        pub fn subtract(self: Self, other: Self) Self {
            var ret = Self.undef();
            for (0..len) |i| {
                ret.data[i] = self.data[i] - other.data[i];
            }
            return ret;
        }

        pub fn negate(self: Self) Self {
            var ret = Self.undef();
            for (0..len) |i| {
                ret.data[i] = -self.data[i];
            }
            return ret;
        }

        pub fn abs2(self: Self) T {
            var ret: T = 0;
            for (0..len) |i| {
                ret += pow(T, self.at(i), 2);
            }
            return ret;
        }

        pub fn scale(self: Self, scalar: T) Self {
            var ret = Self.undef();
            for (0..len) |i| {
                ret.data[i] = self.at(i) * scalar;
            }
            return ret;
        }

        pub fn norm(self: Self) T {
            return sqrt(self.abs2());
        }

        pub fn inner(self: Self, other: Self) T {
            var ret: T = 0;
            for (0..len) |i| {
                ret += self.at(i) * other.at(i);
            }
            return ret;
        }

        pub fn normalize(self: Self) Self {
            return self.scale(1 / self.norm());
        }

        pub fn length(_: Self) usize {
            return len;
        }

        pub fn eq(self: Self, other: Self) bool {
            for (0..len) |i| {
                if (self.at(i) != other.at(i)) {
                    return false;
                }
            }
            return true;
        }

        const IsApproxOptions = struct {
            atol: ?f64 = null,
            rtol: ?f64 = null,
        };

        pub fn isApprox(self: Self, other: Self, opt: IsApproxOptions) bool {
            const atol = opt.atol orelse 0;
            const rtol = opt.rtol orelse sqrt(floatEps(T));
            const n12 = self.subtract(other).norm();
            const n1 = self.norm();
            const n2 = other.norm();
            return n12 <= atol + rtol * @max(n1, n2);
        }
    };
}

test "SVec" {
    const x = SVec(3, f32).fill(1.0);
    const expect = std.testing.expect;
    const expectEqual = std.testing.expectEqual;
    const expectEqualDeep = std.testing.expectEqualDeep;

    try expectEqual(x.at(0), 1.0);
    try expectEqual(x.at(1), 1.0);
    try expectEqual(x.at(2), 1.0);
    try expectEqual(x.length(), 3);

    try expectEqualDeep(x.subtract(x), SVec(3, f32).fill(0));

    try expectEqual(SVec(2, f32).fill(0.0).norm(), 0.0);
    try expectEqual(SVec(2, f32).fill(0.0).abs2(), 0.0);

    try expectEqual(SVec(4, f64).fill(1.0).abs2(), 4.0);
    try expectEqual(SVec(4, f64).fill(2.0).abs2(), 16.0);
    try expectEqual(SVec(4, f64).fill(3.0).norm(), 6.0);

    const y = SVec(3, f32).fill(3.0);
    try expect(!x.eq(y));
    try expect(x.isApprox(x, .{}));
    try expect(!x.isApprox(y, .{}));
    try expect(x.isApprox(y, .{ .atol = 10 }));
    try expect(x.eq(x));
    try expect(x.isApprox(x, .{ .atol = 0, .rtol = 0 }));

    const z = x.add(y);
    try expect(z.eq(SVec(3, f32).fill(4.0)));
    try expect(z.eq(x.scale(4)));

    const a = SVec(3, f32).fromSlice(&[3]f32{ 1.0, 2.0, 3.0 });
    _ = a;
}
