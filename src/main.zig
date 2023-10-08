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
        return (self.start <= x) and (x <= self.stop);
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
        const t2 = res.?[1];
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
            .color = self.color,
        };
    }
};

test "Sphere" {
    const expectEqual = std.testing.expectEqual;
    const expect = std.testing.expect;
    const expectEqualDeep = std.testing.expectEqualDeep;
    {
        const sphere = Sphere{
            .center = Vec.fromXYZ(0, 0, 0),
            .radius = 50,
            .color = RGB.blue,
        };

        const ray = Ray{ .origin = Vec.fromXYZ(0, 0, 10), .velocity = Vec.fromXYZ(0, 0, -5) };
        const res = sphere.hit(ray, Interval{ .start = 0, .stop = std.math.inf(f32) });
        try expect(res != null);
        const hit = res.?;
        try expectEqualDeep(hit, HitRecord{
            .t = 12,
            .point = Vec.fromXYZ(0, 0, -50),
            .normal = Vec.fromXYZ(0, 0, -1),
            .color = RGB.blue,
        });
    }
    {
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
            .color = RGB.red,
        });
    }
}

const HitRecord = struct {
    t: f32,
    point: Vec,
    normal: Vec, // assume normalized
    color: RGB,
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

        var w = std.io.bufferedWriter(file.writer());
        var writer = w.writer();

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
        try w.flush();
    }
};

pub fn convert(comptime T: type, x: anytype) T {
    const ret: T = @floatFromInt(x);
    return ret;
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const alloc = gpa.allocator();
    const img = try Image.initFill(alloc, RGB.white, 1600, 900);
    const dx: f32 = 10 / convert(f32, img.width);
    const dy: f32 = 10 / convert(f32, img.height);

    const hittable = Sphere{
        .center = Vec.fromArray(.{ 0, 0, 0 }),
        .radius = 5,
        .color = RGB.green,
    };
    const z_camera_origin = 10;
    const z_image_plane = 5;
    const ray_origin = Vec.fromArray(.{ 0, 0, z_camera_origin });
    const x_min = -dx * convert(f32, img.width) / 2;
    const y_min = -dy * convert(f32, img.height) / 2;

    // progress
    var i: usize = 0;
    const npixels = img.width * img.height;
    var next_progress_percent: f32 = 0;
    const t_start_gen = std.time.nanoTimestamp();
    for (0..img.height) |iy| {
        for (0..img.width) |ix| {
            i += 1;
            if ((100 * convert(f32, i) / convert(f32, npixels)) >= next_progress_percent) {
                std.debug.print("Progress: {d}%\n", .{next_progress_percent});
                next_progress_percent += 5;
            }
            const x = x_min + convert(f32, ix) * dx;
            const y = y_min + convert(f32, iy) * dy;

            const ray = Ray{
                .origin = ray_origin,
                .velocity = Vec.fromXYZ(x, y, z_image_plane - z_camera_origin),
            };
            const hit = hittable.hit(ray, Interval{ .start = 0, .stop = std.math.inf(f32) }) orelse {
                continue;
            };
            img.at(ix, iy).* = hit.color;
        }
    }
    const seconds_gen = convert(f64, std.time.nanoTimestamp() - t_start_gen) / 1e9;
    std.debug.assert(i == npixels);
    const t_start_save = std.time.nanoTimestamp();
    try img.savePPM("test.ppm");
    const seconds_save = convert(f64, std.time.nanoTimestamp() - t_start_save) / 1e9;
    std.debug.print(
        \\ Done
        \\ seconds generating image : {d}
        \\ seconds saving ppm       : {d}
    , .{ seconds_gen, seconds_save });
}

test {
    std.testing.refAllDecls(svecs);
    std.testing.refAllDecls(@This());
}
