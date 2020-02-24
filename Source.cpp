#include "stdafx.h"
#include <iostream>
#include <vector>
#include <random>
struct Point3D {
	uint16_t x;
	uint16_t y;
	uint16_t z;
	Point3D(uint16_t, uint16_t, uint16_t);
};

Point3D::Point3D(uint16_t x, uint16_t y, uint16_t z) {
	this->x = x;
	this->y = y;
	this->z = z;
}
/* 
struct Vector3D {
	Point3D p;
	Vector3D();
	Vector3D(int, int, int);
	Vector3D(const Vector3D&);
	Vector3D operator+(const Vector3D&);

};

Vector3D::Vector3D() {
	p(0, 0, 0);
}

Vector3D::Vector3D(int x, int y, int z = 0) {
	this->p.x = x;
	this->p.y = y;
	this->p.z = z;
}

Vector3D::Vector3D(const Vector3D& vec) {
	this->p.x = vec.p.x;
	this->p.y = vec.p.y;
	this->p.z = vec.p.z;
}

Vector3D Vector3D::operator+(const Vector3D& vec) {
	return Vector3D(p.x + vec.p.x, p.y + vec.p.y + p.z + vec.p.z);
}
*/
struct Screen {
	int width;
	int height;
	std::vector<Point3D> points;
	Screen();
	Screen(int, int);
	void generatePoints(int);
};
Screen::Screen() {
	width = 800;
	height = 600;
}

Screen::Screen(int w, int h) {
	this->width = w;
	this->height = h;
}

void Screen::generatePoints(int n) {
	points.reserve(n);
	std::random_device rd;
	std::mt19937_64 e2(rd());
	uint16_t upper_bound = 0;
	upper_bound = ~upper_bound;
	std::uniform_int_distribution<uint16_t> res(0, this->width);
	std::uniform_int_distribution<uint16_t> res2(0, this->height);
	for (int i = 0; i < n; i++) {
		uint16_t randomX = res(e2);
		uint16_t randomY = res2(e2);
		Point3D p = Point3D(randomX, randomY, 0);
		points.emplace_back(p);
		// std::cout << random << " ";
	}

}

void radixSort(std::vector<Point3D>& arr) {
	std::vector<std::vector<Point3D> > bucket(2);
	bucket[0] = std::vector<Point3D>(arr.size());
	bucket[1] = std::vector<Point3D>(arr.size());
	int max_length_X = 0;
	int max_length_Y = 0;
	/*
	*	Find largest power of 2 in the array by X and by Y to perform the Radix Sort
	*   o(n)
	*/
	for (int i = 0; i < arr.size(); i++) {
		int log_X = std::log2(arr[i].x);
		int log_Y = std::log2(arr[i].y);
		if (log_X > max_length_X) {
			max_length_X = log_X;
		}
		if (log_Y > max_length_Y) {
			max_length_Y = log_Y;
		}
	}

	/*
	*	Sort by X
	*   O(n*m)
	*/
	for (int i = 0; i < max_length_X; i++) {
		uint16_t mask = 1 << i;
		int idx_bucket = 0;
		for (int j = 0; j < arr.size(); j++) {
			short idx = arr[j].x & mask ? 1 : 0;
			bucket[idx][idx_bucket++] = arr[j];
		}
		idx_bucket = 0;
		for (int c = 0; c < 2; c++) {
			for (int j = 0; j < bucket[c].size(); j++) {
				arr[idx_bucket] = bucket[c][j];
			}
		}
	}
}

int main(){
	Screen* screen = new Screen();
	screen->generatePoints(10);

	delete screen;

}

