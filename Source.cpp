#include <iostream>
#include <vector>
#include <chrono>
#include <cstring>

#include "geometry/Point3D.hpp"
#include "core/Screen.hpp"
#include "core/BVH.hpp"


#define N 250

int main(){
	auto startTime = std::chrono::high_resolution_clock::now();

	Screen* screen = new Screen(800, 600, 600);
	screen->generatePoints(N);
	
	
	
	//for (MortonPrimitive m : mortonPrims) {
	//	std::cout << m.mortonCode << std::endl;
	//}
	
	std::vector<std::shared_ptr<Primitive>> primitives;
	primitives.reserve(N);
	for (int i = 0; i < N; i++) {
		primitives.emplace_back(screen->castPointToPrimitive(i));
	}

	BVH test(primitives, Middle);
	const int boxesToTestCount = 100;
	std::vector<BoundBox> boxesToTest(boxesToTestCount);
	for(int i = 0; i < boxesToTestCount; i++){
        BoundBox boundBox = BoundBox(screen->randomPoint(), screen->randomPoint());
        boxesToTest[i] = boundBox;
	}
	BoundBox boundBox1 = BoundBox(Point3D(0,0,0), Point3D(screen->width, screen->height, screen->depth));
	boxesToTest.emplace_back(boundBox1);
	int yesCount = 0, noCount = 0;
	for(BoundBox& box : boxesToTest){
        if(test.isPrimitiveInBox(box)){
            yesCount++;
        }else{
            noCount++;
        }
	}
	std::cout << "Yes: " << yesCount << " No: " << noCount << std::endl;
	yesCount=0,noCount=0;
	std::cout << "---------------------------------" << std::endl;
	BVH test2(primitives, SAH);
    for(BoundBox& box : boxesToTest){
        if(test2.isPrimitiveInBox(box)){
            yesCount++;
        }else{
            noCount++;
        }
    }

	auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "Yes: " << yesCount << " No: " << noCount << std::endl;
	std::cout << "Time spent: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms\n";
	delete screen;

}

