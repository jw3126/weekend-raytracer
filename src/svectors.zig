const std = @import("std");
const sqrt = std.math.sqrt;
const floatEps = std.math.floatEps;
const pow = std.math.pow;

pub fn Vec(comptime len: comptime_int, comptime T: type) type {
    return struct {
        const Self = @This();
        data: [len]T,

        fn undef() Self {
            var ret = Self{
                .data = undefined,
            };
            return ret;
        }

        pub fn unitvector(i: usize) Self {
            var ret = Self.fill(0);
            ret.data[i] = 1;
            return ret;
        }

        pub fn fromXYZ(x: T, y: T, z: T) Self {
            comptime if (len != 3) {
                @compileError("fromXYZ requires 3D vectors");
            };
            var ret = Self.undef();
            ret.data[0] = x;
            ret.data[1] = y;
            ret.data[2] = z;
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

        pub fn fromArray(array: [len]T) Self {
            return Self.fromSlice(&array);
        }

        pub fn at(self: Self, index: anytype) T {
            return self.data[index];
        }

        pub fn map2(self: Self, other: Self, comptime f: fn (T, T) T) Self {
            var ret = Self.undef();
            for (0..len) |i| {
                ret.data[i] = f(self.at(i), other.at(i));
            }
            return ret;
        }

        pub fn map(self: Self, comptime f: fn (T) T) Self {
            var ret = Self.undef();
            for (0..len) |i| {
                ret.data[i] = f(self.at(i));
            }
            return ret;
        }

        fn add_scalar(x: T, y: T) T {
            return x + y;
        }

        fn sub_scalar(x: T, y: T) T {
            return x - y;
        }

        fn negate_scalar(x: T) T {
            return -x;
        }

        pub fn add(self: Self, other: Self) Self {
            return self.map2(other, add_scalar);
        }

        pub fn subtract(self: Self, other: Self) Self {
            return self.map2(other, sub_scalar);
        }

        pub fn negate(self: Self) Self {
            return self.map(negate_scalar);
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

        pub fn cross(self: Self, other: Self) Self {
            comptime if (len != 3) {
                @compileError("cross product requires 3D vectors");
            };
            return Self{
                .data = [3]T{
                    self.at(1) * other.at(2) - self.at(2) * other.at(1),
                    self.at(2) * other.at(0) - self.at(0) * other.at(2),
                    self.at(0) * other.at(1) - self.at(1) * other.at(0),
                },
            };
        }
    };
}

test "Vec" {
    const x = Vec(3, f32).fill(1.0);
    const expect = std.testing.expect;
    const expectEqual = std.testing.expectEqual;
    const expectEqualDeep = std.testing.expectEqualDeep;

    try expectEqual(x.at(0), 1.0);
    try expectEqual(x.at(1), 1.0);
    try expectEqual(x.at(2), 1.0);
    try expectEqual(x.length(), 3);

    try expectEqualDeep(x.subtract(x), Vec(3, f32).fill(0));

    try expectEqual(Vec(2, f32).fill(0.0).norm(), 0.0);
    try expectEqual(Vec(2, f32).fill(0.0).abs2(), 0.0);

    try expectEqual(Vec(4, f64).fill(1.0).abs2(), 4.0);
    try expectEqual(Vec(4, f64).fill(2.0).abs2(), 16.0);
    try expectEqual(Vec(4, f64).fill(3.0).norm(), 6.0);

    const y = Vec(3, f32).fill(3.0);
    try expect(!x.eq(y));
    try expect(x.isApprox(x, .{}));
    try expect(!x.isApprox(y, .{}));
    try expect(x.isApprox(y, .{ .atol = 10 }));
    try expect(x.eq(x));
    try expect(x.isApprox(x, .{ .atol = 0, .rtol = 0 }));

    const z = x.add(y);
    try expect(z.eq(Vec(3, f32).fill(4.0)));
    try expect(z.eq(x.scale(4)));

    const a = Vec(3, f32).fromArray(.{ 1.0, 2.0, 3.0 });
    try expectEqual(a.at(0), 1.0);
    try expectEqual(a.at(1), 2.0);
    try expectEqual(a.at(2), 3.0);
}

test "cross" {
    const expectEqualDeep = std.testing.expectEqualDeep;
    const e0 = Vec(3, u8).unitvector(0);
    const e1 = Vec(3, u8).unitvector(1);
    const e2 = Vec(3, u8).unitvector(2);
    try expectEqualDeep(e0.cross(e1), e2);
    try expectEqualDeep(e1.cross(e2), e0);
    try expectEqualDeep(e2.cross(e0), e1);
}
