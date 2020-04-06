#include "Point3D.hpp"
#pragma once

struct Vector3D {
	int x;
	int y;
	int z;

	Vector3D(Point3D&);
	Vector3D();
	Vector3D diagonal();
};

Vector3D::Vector3D(Point3D& p) {
	x = p.x;
	y = p.y;
	z = p.z;
}

