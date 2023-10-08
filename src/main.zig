const std = @import("std");
const Allocator = std.mem.Allocator;
const svecs = @import("svectors.zig");
const Vec = svecs.Vec(3, f32);

const Ray = struct {
    origin: Vec,
    velocity: Vec,

    pub fn at(self: Ray, t: f32) Vec {
        return self.origin.add(self.velocity.scale(t));
    }
};

pub fn abc_formula(a: f32, b: f32, c: f32) ?[2]f32 {
    if (a == 0) {
        // should we return a single solution?
        // should we error?
        return null;
    }
    const disc2 = b * b - 4 * a * c;
    if (disc2 < 0) {
        return null;
    }
    const disc = @sqrt(disc2);
    const sol1 = (-b - disc) / (2 * a);
    const sol2 = (-b + disc) / (2 * a);
    if (a > 0) {
        return .{ sol1, sol2 };
    } else {
        return .{ sol2, sol1 };
    }
}

test "abc forumla" {
    const expectEqual = std.testing.expectEqual;
    try expectEqual(abc_formula(1, 0, 0), .{ 0, 0 });
    try expectEqual(abc_formula(0, 1, -1), null);
    try expectEqual(abc_formula(1, 0, 1), null);
    try expectEqual(abc_formula(1, 0, -1), .{ -1, 1 });
    try expectEqual(abc_formula(1, 0, -1), .{ -1, 1 });
    try expectEqual(abc_formula(1, -3, 2), .{ 1, 2 });
    try expectEqual(abc_formula(2, -6, 4), .{ 1, 2 });
    try expectEqual(abc_formula(-1, 3, -2), .{ 1, 2 });
}

const Interval = struct {
    const Self = @This();
    start: f32,
    stop: f32,
    pub fn contains(self: Self, x: f32) bool {
        return x >= self.start and x <= self.stop;
    }
};

const Sphere = struct {
    center: Vec,
    radius: f32,
    color: RGB,

    pub fn hit(self: Sphere, ray: Ray, ts: Interval) ?HitRecord {
        const oc = ray.origin.subtract(self.center);
        const v = ray.velocity;
        const a = v.abs2();
        const b = 2 * oc.inner(v);
        const c = oc.inner(oc) - self.radius * self.radius;
        const res = abc_formula(a, b, c);
        if (res == null) {
            return null;
        }
        const t1 = res.?[0];
        const t2 = res.?[0];
        var t = std.math.nan(f32);
        if (ts.contains(t1)) {
            t = t1;
        } else if (ts.contains(t2)) {
            t = t2;
        } else {
            return null;
        }
        const point = ray.at(t);
        const normal = point.subtract(self.center).scale(1 / self.radius);
        return HitRecord{
            .t = t,
            .point = point,
            .normal = normal,
        };
    }
};

test "Sphere" {
    const expectEqual = std.testing.expectEqual;
    const expect = std.testing.expect;
    const expectEqualDeep = std.testing.expectEqualDeep;
    const sphere = Sphere{
        .center = Vec.fromArray(.{ 0, 0, 0 }),
        .radius = 1,
        .color = RGB.red,
    };
    const ray = Ray{
        .origin = Vec.fromArray(.{ 0, 0, 10 }),
        .velocity = Vec.fromArray(.{ 0, 0, -1 }),
    };
    const res = sphere.hit(ray, Interval{ .start = 0, .stop = 10 });
    try expect(res != null);
    const hit = res.?;
    try expectEqual(hit.t, 9);
    try expectEqualDeep(hit, HitRecord{
        .t = 9,
        .point = Vec.fromArray(.{ 0, 0, 1 }),
        .normal = Vec.fromArray(.{ 0, 0, 1 }),
    });
}

const HitRecord = struct {
    t: f32,
    point: Vec,
    normal: Vec, // assume normalized
    // material: *const Sphere,
};

const RGB = struct {
    const Self = @This();
    r: f32,
    g: f32,
    b: f32,
    const green = Self{ .r = 0, .g = 1, .b = 0 };
    const blue = Self{ .r = 0, .g = 0, .b = 1 };
    const red = Self{ .r = 1, .g = 0, .b = 0 };
    const black = Self{ .r = 0, .g = 0, .b = 0 };
    const white = Self{ .r = 1, .g = 1, .b = 1 };
};

const Image = struct {
    const Self = @This();
    width: usize,
    height: usize,
    pixelValues: []RGB, // SOA?
    alloc: Allocator,

    pub fn initFill(alloc: Allocator, color: RGB, width: usize, height: usize) !Self {
        var pixelValues = try alloc.alloc(RGB, width * height);
        for (pixelValues) |*p| {
            p.* = color;
        }
        return Self{
            .width = width,
            .height = height,
            .pixelValues = pixelValues,
            .alloc = alloc,
        };
    }

    pub fn at(self: Self, x: usize, y: usize) *RGB {
        const i = self.toLinearIndex(x, y);
        return &self.pixelValues[i];
    }

    pub fn setAt(self: *Self, x: usize, y: usize, color: RGB) void {
        const i = self.toLinearIndex(x, y);
        self.pixelValues[i] = color;
    }

    fn toLinearIndex(self: Self, x: usize, y: usize) usize {
        return y * self.width + x;
    }

    pub fn deinit(self: *Self) void {
        self.alloc.free(self.pixelValues);
    }

    pub fn savePPM(self: Self, path: []const u8) !void {
        const file = try std.fs.cwd().createFile(path, .{});
        defer file.close();
        const writer = file.writer();
        try writer.print("P3\n", .{});
        try writer.print("{d} {d}\n", .{ self.width, self.height });
        try writer.print("255\n", .{}); // max color, do we want to adjust?
        for (0..self.height) |iy| {
            for (0..self.width) |ix| {
                const color = self.at(ix, iy);
                const r: u8 = @intFromFloat(std.math.round(color.r * 255));
                const g: u8 = @intFromFloat(std.math.round(color.g * 255));
                const b: u8 = @intFromFloat(std.math.round(color.b * 255));
                try writer.print("{d} {d} {d}\n", .{ r, g, b });
            }
        }
    }
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const alloc = gpa.allocator();
    const img = try Image.initFill(alloc, RGB.white, 256, 256);
    for (0..img.height) |y| {
        for (0..img.width) |x| {
            const fx: f32 = @floatFromInt(x);
            const fy: f32 = @floatFromInt(y);
            const fz: f32 = @floatFromInt(try std.math.mod(usize, x * y, 255));
            img.at(x, y).* = RGB{ .r = fx / 255, .g = fy / 255, .b = fz / 255 };
        }
    }
    try img.savePPM("test.ppm");
}

test {
    std.testing.refAllDecls(svecs);
    std.testing.refAllDecls(@This());
}
