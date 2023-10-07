const std = @import("std");
const Allocator = std.mem.Allocator;
const svecs = @import("svectors.zig");
const SVec = svecs.SVec;

const RGB = struct {
    const Self = @This();
    r: f32,
    g: f32,
    b: f32,
    const green = Self{ .r = 0, .g = 1, .b = 0 };
    const blue = Self{ .r = 0, .g = 0, .b = 1 };
    const red = Self{ .r = 1, .g = 0, .b = 0 };
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
    const img = try Image.initFill(alloc, RGB.blue, 256, 256);
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
