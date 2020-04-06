#pragma once
#include <vector>
#include <random>
#include "Point3D.hpp"
#include "BoundBox.hpp"


struct Primitive {
	Primitive(Point3D p) {
		box = BoundBox(p - 10, p + 10);
	}
	Primitive() {}
	BoundBox box;
};

struct Screen {
	int width;
	int height;
	int depth;
	std::vector<Point3D> points;
	Screen();
	Screen(int, int, int);
	void generatePoints(int);
	void showPoints();
	Primitive* castPointToPrimitive(int);
};
Screen::Screen() {
	width = 800;
	height = 600;
	depth = 300;
}

Screen::Screen(int w, int h, int d) {
	this->width = w;
	this->height = h;
	this->depth = d;
}

void Screen::generatePoints(int n) {
	points.reserve(n);
	std::random_device rd;
	std::mt19937_64 e2(rd());
	uint32_t upper_bound = 0;
	upper_bound = ~upper_bound;
	std::uniform_int_distribution<uint32_t> res(1, this->width-1);
	std::uniform_int_distribution<uint32_t> res2(1, this->height-1);
	std::uniform_int_distribution<uint32_t> res3(1, this->depth-1);
	for (int i = 0; i < n; i++) {
		uint32_t randomX = res(e2);
		uint32_t randomY = res2(e2);
		uint32_t randomZ = res3(e2);
		Point3D p = Point3D(randomX, randomY, randomZ);
		points.emplace_back(p);
		// std::cout << random << " ";
	}

}

void Screen::showPoints() {
	for (Point3D p : points) {
		printf("\tX: %f, Y: %f, Z: %f\n", p.x, p.y, p.z);
	}
}

Primitive* Screen::castPointToPrimitive(int index) {
	Primitive* primitive = new Primitive(this->points[index]);
	return primitive;
}

