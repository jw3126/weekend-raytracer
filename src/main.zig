const std = @import("std");
const Allocator = std.mem.Allocator;
const svecs = @import("svectors.zig");
const Vec = svecs.Vec(3, f64);

const Ray = struct {
    origin: Vec,
    velocity: Vec,

    pub fn at(self: Ray, t: f64) Vec {
        return self.origin.add(self.velocity.scale(t));
    }
};

const green = Vec.fromXYZ(0, 1, 0);
const blue = Vec.fromXYZ(0, 0, 1);
const red = Vec.fromXYZ(1, 0, 0);
const black = Vec.fromXYZ(0, 0, 0);
const white = Vec.fromXYZ(1, 1, 1);
const grey = Vec.fromXYZ(0.5, 0.5, 0.5);

pub fn abc_formula(a: anytype, b: anytype, c: anytype) ?[2]@TypeOf(a, b, c) {
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

fn abc_formula64(a: f64, b: f64, c: f64) ?[2]f64 {
    return abc_formula(a, b, c);
}

fn abc_formula32(a: f32, b: f32, c: f32) ?[2]f32 {
    return abc_formula(a, b, c);
}

test "abc forumla" {
    const expectEqual = std.testing.expectEqual;

    try expectEqual(abc_formula64(1, 0, 0), .{ 0, 0 });
    try expectEqual(abc_formula64(0, 1, -1), null);
    try expectEqual(abc_formula64(1, 0, 1), null);
    try expectEqual(abc_formula64(1, 0, -1), .{ -1, 1 });
    try expectEqual(abc_formula64(1, 0, -1), .{ -1, 1 });
    try expectEqual(abc_formula32(1, -3, 2), .{ 1, 2 });
    try expectEqual(abc_formula32(2, -6, 4), .{ 1, 2 });
    try expectEqual(abc_formula32(-1, 3, -2), .{ 1, 2 });
}

const Interval = struct {
    const Self = @This();
    start: f64,
    stop: f64,
    pub fn contains(self: Self, x: f64) bool {
        return (self.start <= x) and (x <= self.stop);
    }

    pub fn clamp(self: Self, x: f64) f64 {
        return std.math.clamp(x, self.start, self.stop);
    }
};

const Sphere = struct {
    center: Vec,
    radius: f64,
    color: Vec,

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
        var t = std.math.nan(f64);
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
            .color = blue,
        };

        const ray = Ray{ .origin = Vec.fromXYZ(0, 0, 10), .velocity = Vec.fromXYZ(0, 0, -5) };
        const res = sphere.hit(ray, Interval{ .start = 0, .stop = std.math.inf(f64) });
        try expect(res != null);
        const hit = res.?;
        try expectEqualDeep(hit, HitRecord{
            .t = 12,
            .point = Vec.fromXYZ(0, 0, -50),
            .normal = Vec.fromXYZ(0, 0, -1),
            .color = blue,
        });
    }
    {
        const sphere = Sphere{
            .center = Vec.fromArray(.{ 0, 0, 0 }),
            .radius = 1,
            .color = red,
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
            .color = red,
        });
    }
}

const HitRecord = struct {
    t: f64,
    point: Vec,
    normal: Vec, // assume normalized
    color: Vec,
    // material: *const Sphere,
};

const Image = struct {
    const Self = @This();
    width: usize,
    height: usize,
    pixelValues: []Vec, // SOA?
    alloc: Allocator,

    pub fn initFill(alloc: Allocator, color: Vec, width: usize, height: usize) !Self {
        var pixelValues = try alloc.alloc(Vec, width * height);
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

    pub fn at(self: Self, x: usize, y: usize) *Vec {
        const i = self.toLinearIndex(x, y);
        return &self.pixelValues[i];
    }

    pub fn setAt(self: *Self, x: usize, y: usize, color: Vec) void {
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
                const r: u8 = Self.ppmScalar(color.at(0));
                const g: u8 = Self.ppmScalar(color.at(1));
                const b: u8 = Self.ppmScalar(color.at(2));
                try writer.print("{d} {d} {d}\n", .{ r, g, b });
            }
        }
        try w.flush();
    }

    fn gammaCorrect(x: f64) f64 {
        return std.math.sqrt(x);
    }

    fn ppmScalar(x: f64) u8 {
        const g = gammaCorrect(x);
        const y: u8 = @intFromFloat(std.math.round(g * 255));
        return y;
    }
};

const RayTraceable = struct {
    const Self = @This();
    hittables: []const Sphere,
    max_reflections: u32 = 1,

    fn calcFirstHit(self: RayTraceable, ray: Ray, ts: Interval) ?HitRecord {
        var start: f64 = ts.start;
        var stop: f64 = ts.stop;
        var maybe_hit: ?HitRecord = null;
        for (self.hittables) |*hittable| {
            if (hittable.hit(ray, .{ .start = start, .stop = stop })) |hit| {
                stop = hit.t;
                maybe_hit = hit;
            }
        }
        return maybe_hit;
    }

    fn rayColor(
        self: RayTraceable,
        rng: std.rand.Random,
        ray: Ray,
    ) Vec {
        var ret = Vec.fromXYZ(0, 0, 0);
        var weight: f64 = 1.0;
        var r = ray;
        for (0..self.max_reflections) |i| {
            _ = i;
            const start: f64 = 1e-3;
            const stop: f64 = std.math.inf(f64);
            const maybe_hit = self.calcFirstHit(r, .{ .start = start, .stop = stop });
            if (maybe_hit == null) {
                ret = ret.add(backgroundColor(r).scale(weight));
            } else {
                ret = ret.add((maybe_hit.?.color).scale(1.0 * weight));
                weight *= 0.0;
                const new_direction = Self.randomOnHemisphere(rng, r.velocity);
                const new_origin = maybe_hit.?.point;
                r = Ray{ .origin = new_origin, .velocity = new_direction };
            }
        }
        return ret;
    }

    fn randomOnHemisphere(rng: std.rand.Random, v: Vec) Vec {
        while (true) {
            const x = rng.float(f64) * 2.0 - 1.0;
            const y = rng.float(f64) * 2.0 - 1.0;
            const z = rng.float(f64) * 2.0 - 1.0;
            const candidate = Vec.fromXYZ(x, y, z);
            if (v.inner(candidate) > 0) {
                return candidate.normalize();
            }
        }
    }

    fn backgroundColor(ray: Ray) Vec {
        const y = 0.5 * (ray.velocity.at(1) + 1.0);
        return Vec.lerp(white, Vec.fromXYZ(0.5, 0.7, 1.0), y);
    }
};

pub fn convert(comptime Target: type, x: anytype) Target {
    const ret: Target = @floatFromInt(x);
    return ret;
}

pub fn uniform(comptime T: type, rng: std.rand.Random, min: anytype, max: anytype) T {
    return rng.float(T) * (max - min) + min;
}

pub fn main() !void {
    const T = f64;
    const t_start = std.time.nanoTimestamp();
    var prng = std.rand.DefaultPrng.init(42);
    const rng: std.rand.Random = prng.random();

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const alloc = gpa.allocator();
    const img = try Image.initFill(alloc, white, 1600 / 2, 900 / 2);
    const aspect_ratio = convert(T, img.width) / convert(T, img.height);

    const s1 = Sphere{
        .center = Vec.fromArray(.{ 0, 0, -1 }),
        .radius = 0.5,
        .color = red,
    };

    const s_big = Sphere{
        .center = Vec.fromArray(.{ 0, -1000, 100 }),
        .radius = 10,
        .color = grey,
    };

    const hittables: []const Sphere = &[_]Sphere{ s1, s_big };

    const traceable = RayTraceable{
        .hittables = hittables,
    };
    const camera_origin = Vec.fromXYZ(0, 0, 0);
    const focal_length = 1;
    const viewport_height = 2;
    const viewport_width = viewport_height * aspect_ratio;
    const viewport_u = Vec.fromXYZ(viewport_width, 0, 0);
    const viewport_v = Vec.fromXYZ(0, -viewport_height, 0);

    const pixel_delta_u = viewport_u.scale(1 / convert(T, img.width));
    const pixel_delta_v = viewport_v.scale(1 / convert(T, img.height));

    const viewport_upper_left = camera_origin
        .add(Vec.fromXYZ(0, 0, -focal_length))
        .add(viewport_u.scale(-0.5))
        .add(viewport_v.scale(-0.5));

    const pixel00_loc = viewport_upper_left
        .add(pixel_delta_u.scale(0.5))
        .add(pixel_delta_v.scale(0.5));
    // progress
    var i: usize = 0;
    const npixels = img.width * img.height;
    var next_progress_percent: T = 0;
    const t_start_gen = std.time.nanoTimestamp();
    const nsamples = 10;
    for (0..img.height) |iv| {
        for (0..img.width) |iu| {
            i += 1;
            if ((100 * convert(T, i) / convert(T, npixels)) >= next_progress_percent) {
                std.debug.print("Progress: {d}%\n", .{next_progress_percent});
                next_progress_percent += 5;
            }

            const pixel_center = pixel00_loc
                .add(pixel_delta_u.scale(convert(T, iu)))
                .add(pixel_delta_v.scale(convert(T, iv)));

            var color = Vec.fill(0);
            for (0..nsamples) |_| {
                const pixel_sample = pixel_center
                    .add(pixel_delta_u.scale(uniform(T, rng, -0.5, 0.5)))
                    .add(pixel_delta_v.scale(uniform(T, rng, -0.5, 0.5)));
                const ray = Ray{
                    .origin = camera_origin,
                    .velocity = pixel_sample.subtract(camera_origin).normalize(),
                };
                const hit_color = traceable.rayColor(rng, ray);
                color = color.add(hit_color);
            }
            img.at(iu, iv).* = color.scale(convert(T, 1) / nsamples);
        }
    }
    const seconds_gen = convert(f64, std.time.nanoTimestamp() - t_start_gen) / 1e9;
    std.debug.assert(i == npixels);
    const t_start_save = std.time.nanoTimestamp();
    try img.savePPM("test.ppm");
    const seconds_save = convert(f64, std.time.nanoTimestamp() - t_start_save) / 1e9;
    const seconds_total = convert(f64, std.time.nanoTimestamp() - t_start) / 1e9;
    std.debug.print(
        \\ Done
        \\ seconds generating image : {d}
        \\ seconds saving ppm       : {d}
        \\ seconds total            : {d}
    , .{ seconds_gen, seconds_save, seconds_total });
}

test {
    std.testing.refAllDecls(svecs);
    std.testing.refAllDecls(@This());
}
