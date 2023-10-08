all:
	zig build run
	pnm2png test.ppm > test.png
	xdg-open test.png
