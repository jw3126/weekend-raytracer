const std = @import("std");
const Allocator = std.mem.Allocator;
const svecs = @import("svectors.zig");
const SVec = svecs.SVec;

const RGB = struct { r: f32, g: f32, b: f32 };

const Image = struct {
    const Self = @This();
    width: u32,
    height: u32,
    pixelValues: []RGB, // SOA?
    alloc: Allocator,

    pub fn initFill(alloc: Allocator, color: RGB, width: u32, height: u32) Self {
        _ = color;
        var pixelValues = try alloc.alloc(f32, width * height);
        return Self{
            .width = width,
            .height = height,
            .pixelValues = pixelValues,
            .alloc = alloc,
        };
    }

    pub fn at(self: Self, x: u32, y: u32) RGB {
        const i = self.toLinearIndex(x, y);
        return self.pixelValues[i];
    }

    pub fn setAt(self: *Self, x: u32, y: u32, color: RGB) void {
        const i = self.toLinearIndex(x, y);
        self.pixelValues[i] = color;
    }

    fn toLinearIndex(self: Self, x: u32, y: u32) u32 {
        return y * self.width + x;
    }

    pub fn deinit(self: *Self) void {
        self.alloc.free(self.pixelValues);
    }

    pub fn savePPM(path: []const u8, self: Self) !void {
        const file = try std.fs.cwd().createFile(path, .{});
        defer file.close();
        try file.writeAll("P3\n");
        try file.writeAll(&[_]u8{ self.width, ' ', self.height });
    }
};

pub fn main() !void {
    // Prints to stderr (it's a shortcut based on `std.io.getStdErr()`)
    std.debug.print("All your {s} are belong to us.\n", .{"codebase"});

    // stdout is for the actual output of your application, for example if you
    // are implementing gzip, then only the compressed bytes should be sent to
    // stdout, not any debugging messages.
    const stdout_file = std.io.getStdOut().writer();
    var bw = std.io.bufferedWriter(stdout_file);
    const stdout = bw.writer();

    try stdout.print("Run `zig build test` to run the tests.\n", .{});

    try bw.flush(); // don't forget to flush!
}

test {
    std.testing.refAllDecls(svecs);
    std.testing.refAllDecls(@This());
}
