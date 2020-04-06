#include <iostream>
#include <vector>
#include <chrono>

#include "Point3D.hpp"
#include "Screen.hpp"
#include "BVH.hpp"


#define N 150

int main(){
	auto startTime = std::chrono::high_resolution_clock::now();

	Screen* screen = new Screen(800, 600, 300);
	screen->generatePoints(N);
	
	
	
	//for (MortonPrimitive m : mortonPrims) {
	//	std::cout << m.mortonCode << std::endl;
	//}
	
	std::vector<std::shared_ptr<Primitive>> primitives;
	primitives.reserve(N);
	for (int i = 0; i < N; i++) {
		primitives.emplace_back(screen->castPointToPrimitive(i));
	}

	BVH test(primitives);
	auto endTime = std::chrono::high_resolution_clock::now();

	std::cout << "Time spent: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms\n";

	getchar();
	delete screen;

}

