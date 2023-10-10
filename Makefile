all:
	zig build run -Doptimize=ReleaseFast
	pnm2png test.ppm > test.png
	xdg-open test.png
