#pragma once
#include <cstdint>

struct Point3D {
	float x;
	float y;
	float z;
	Point3D(float, float, float);
	Point3D();

	float operator[](int);
	float operator[](int) const;
	Point3D operator+(int);
	Point3D operator-(int);
	Point3D operator-(Point3D&);
    bool operator==(const Point3D &p) const;
};

Point3D::Point3D() {
	x = 0;
	y = 0;
	z = 0;
}

Point3D::Point3D(float x, float y, float z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

bool operator<(Point3D a, Point3D b) {
	float xSquare = a.x * a.x;
	float ySquare = a.y * a.y;
	float zSquare = a.z * a.z;

	float x2Square = b.x * b.x;
	float y2Square = b.y * b.y;
	float z2Square = b.z * b.z;

	float sum = std::sqrt(xSquare + ySquare + z2Square) - std::sqrt(x2Square + y2Square + z2Square);
	return sum < 0 ||
		sum == 0 && xSquare < x2Square ||
		sum == 0 && xSquare == x2Square && ySquare < y2Square ||
		sum == 0 && xSquare == x2Square && ySquare == y2Square && zSquare < z2Square;
}

bool operator>(Point3D a, Point3D b) {
	float xSquare = a.x * a.x;
	float ySquare = a.y * a.y;
	float zSquare = a.z * a.z;

	float x2Square = b.x * b.x;
	float y2Square = b.y * b.y;
	float z2Square = b.z * b.z;

	float sum = std::sqrt(xSquare + ySquare + z2Square) - std::sqrt(x2Square + y2Square + z2Square);
	return sum > 0 ||
		sum == 0 && xSquare > x2Square ||
		sum == 0 && xSquare == x2Square && ySquare > y2Square ||
		sum == 0 && xSquare == x2Square && ySquare == y2Square && zSquare > z2Square;
}

float Point3D::operator[](const int i) {
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

float Point3D::operator[](const int i) const {
	if (i == 0) return x;
	if (i == 1) return y;
	return z;
}

bool Point3D::operator==(const Point3D &p) const {
    return x == p.x && y == p.y && z == p.z;
}
