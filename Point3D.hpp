#include <cstdint>
#pragma once
struct Point3D {
	float x;
	float y;
	float z;
	Point3D(uint32_t, uint32_t, uint32_t);
	Point3D();

	int operator[](int);
	int operator[](int) const;
	Point3D operator+(int);
	Point3D operator-(int);
	Point3D operator-(Point3D&);
};

Point3D::Point3D() {
	x = 0;
	y = 0;
	z = 0;
}

Point3D::Point3D(uint32_t x, uint32_t y, uint32_t z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

bool operator<(Point3D a, Point3D b) {
	uint32_t xSquare = a.x * a.x;
	uint32_t ySquare = a.y * a.y;
	uint32_t zSquare = a.z * a.z;

	uint32_t x2Square = b.x * b.x;
	uint32_t y2Square = b.y * b.y;
	uint32_t z2Square = b.z * b.z;

	int64_t sum = std::sqrt(xSquare + ySquare + z2Square) - std::sqrt(x2Square + y2Square + z2Square);
	return sum < 0 ||
		sum == 0 && xSquare < x2Square ||
		sum == 0 && xSquare == x2Square && ySquare < y2Square ||
		sum == 0 && xSquare == x2Square && ySquare == y2Square && zSquare < z2Square;
}

bool operator>(Point3D a, Point3D b) {
	uint32_t xSquare = a.x * a.x;
	uint32_t ySquare = a.y * a.y;
	uint32_t zSquare = a.z * a.z;

	uint32_t x2Square = b.x * b.x;
	uint32_t y2Square = b.y * b.y;
	uint32_t z2Square = b.z * b.z;

	int32_t sum = std::sqrt(xSquare + ySquare + z2Square) - std::sqrt(x2Square + y2Square + z2Square);
	return sum > 0 ||
		sum == 0 && xSquare > x2Square ||
		sum == 0 && xSquare == x2Square && ySquare > y2Square ||
		sum == 0 && xSquare == x2Square && ySquare == y2Square && zSquare > z2Square;
}

int Point3D::operator[](int i) {
	if (i == 0) return x;
	if (i == 1) return y;
	return z;
}

Point3D Point3D::operator+(int i) {
	this->x += i;
	this->y += i;
	this->z += i;
	return *this;
}

Point3D Point3D::operator-(const int i) {
	this->x -= i;
	this->y -= i;
	this->z -= i;
	return *this;
}

Point3D Point3D::operator-(Point3D& p) {
	this->x -= p.x;
	this->y -= p.y;
	this->z -= p.z;
	return *this;
}

int Point3D::operator[](const int i) const {
	if (i == 0) return x;
	if (i == 1) return y;
	return z;
}

