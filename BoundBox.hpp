#include "Point3D.hpp"
#include "Vector3D.hpp"
#pragma once


struct BoundBox { 
	Point3D pMin;
	Point3D pMax;

	BoundBox(Point3D);
	BoundBox(Point3D, Point3D);
	BoundBox();
	void setBounds(BoundBox);
	void Union(BoundBox);

	BoundBox Union(BoundBox&, Point3D&);
	BoundBox Union(BoundBox, BoundBox);
	BoundBox unite(BoundBox, BoundBox);
	BoundBox unite(BoundBox);

	const Point3D offset(const Point3D&);
	Point3D diagonal();
	const int MaximumExtent();
	float surfaceArea();
};


BoundBox::BoundBox() {
	float minNum = 0;
	pMin = Point3D(800, 600, 300);
	pMax = Point3D(minNum, minNum, minNum);
}
BoundBox::BoundBox(Point3D p){
	pMin = p;
	pMax = p;
}

BoundBox::BoundBox(Point3D p1, Point3D p2) {
	pMin = Point3D(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
	pMax = Point3D(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
		
}

BoundBox BoundBox::Union(BoundBox& box, Point3D& p) {
	BoundBox newBox;
	newBox.pMin = Point3D(std::min(box.pMin.x, p.x), std::min(box.pMin.y, p.y), std::min(box.pMin.z, p.z));
	newBox.pMax = Point3D(std::max(box.pMax.x, p.x), std::max(box.pMax.y, p.y), std::max(box.pMax.z, p.z));
	return newBox;
}

BoundBox BoundBox::Union(BoundBox box1, BoundBox box2) {
	BoundBox newBox;
	newBox.pMin = std::min(box1.pMin, box2.pMin);
	newBox.pMax = std::max(box1.pMax, box2.pMax);
	return newBox;
}

BoundBox Union(BoundBox box1, BoundBox box2) {
	BoundBox newBox;
	newBox.pMin = std::min(box1.pMin, box2.pMin);
	newBox.pMax = std::max(box1.pMax, box2.pMax);
	return newBox;
}

BoundBox BoundBox::unite(BoundBox b1, BoundBox b2) {
	bool x = (b1.pMax.x >= b2.pMin.x) && (b1.pMin.x <= b2.pMax.x);
	bool y = (b1.pMax.y >= b2.pMin.y) && (b1.pMin.y <= b2.pMax.y);
	bool z = (b1.pMax.z >= b2.pMin.z) && (b1.pMin.z <= b2.pMax.z);
	if (x && y && z) {
		return Union(b1, b2);
	}
}

BoundBox BoundBox::unite(BoundBox b2) {
	bool x = (this->pMax.x >= b2.pMin.x) && (this->pMin.x <= b2.pMax.x);
	bool y = (this->pMax.y >= b2.pMin.y) && (this->pMin.y <= b2.pMax.y);
	bool z = (this->pMax.z >= b2.pMin.z) && (this->pMin.z <= b2.pMax.z);
	if (x && y && z) {
		return Union(*this, b2);
	}
	else return *this;
}

const int BoundBox::MaximumExtent() {
	Point3D d = Point3D(this->pMax.x - this->pMin.x, this->pMax.y - this->pMin.y, this->pMax.z - this->pMin.z); // diagonal
	if (d.x > d.y && d.x > d.z) {
		return 0;
	}
	else if (d.y > d.z) {
		return 1;
	}
	else {
		return 2;
	}
}

float BoundBox::surfaceArea() {
	Point3D d = Point3D(this->pMax.x - this->pMin.x, this->pMax.y - this->pMin.y, this->pMax.z - this->pMin.z); // diagonal
	return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
}

const Point3D BoundBox::offset(const Point3D& p) {
	Point3D o = Point3D(p.x - pMin.x, p.y - pMin.y, p.z - pMin.z);
	
	if (pMax.x > pMin.x) o.x /= pMax.x - pMin.x;
	if (pMax.y > pMin.y) o.y /= pMax.y - pMin.y;
	if (pMax.z > pMin.z) o.z /= pMax.z - pMin.z;
	return o;
}