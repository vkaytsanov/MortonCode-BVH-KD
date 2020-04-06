#pragma once
#include <vector>
#include <cassert>
#include <algorithm>
#include "memory.hpp"
#include "Screen.hpp"
#include "Point3D.hpp"
#include "BoundBox.hpp"


#pragma once
enum Axis{
	X, Y, Z
};
struct MortonPrimitive{
	int primitiveIndex;
	uint32_t mortonCode;
};

struct BVHPrimitiveInfo {
	BVHPrimitiveInfo() {}
	BVHPrimitiveInfo(int primitiveNumber, const BoundBox& box) : primitiveNumber(primitiveNumber), box(box),
		centroid(Point3D(box.pMin.x* 0.5f + box.pMax.x * 0.5f, box.pMin.y* 0.5f + box.pMax.y * 0.5f, box.pMin.z* 0.5f + box.pMax.z * 0.5f)) {}

	int primitiveNumber;
	BoundBox box;
	Point3D centroid;
};

struct BVHNode {
	void InitLeaf(int first, int n, const BoundBox& b) {
		firstPrimOffset = first;
		nPrimitives = n;
		box = b;
		children[0] = children[1] = nullptr;

	}
	void InitInterior(int axis, BVHNode* c0, BVHNode* c1) {
		//assert(c0 != nullptr || c1 != nullptr);
		children[0] = c0;
		children[1] = c1;
		box = Union(c0->box, c1->box);
		splitAxis = axis;
		nPrimitives = 0;
	}
	BoundBox box;
	BVHNode* children[2];
	int splitAxis, firstPrimOffset, nPrimitives;
};

struct LinearBVHNode {
	BoundBox bounds;
	union {
		int primitivesOffset;   // leaf
		int secondChildOffset;  // interior
	};
	uint16_t nPrimitives;  // 0 -> interior node
	uint8_t axis;          // interior node: xyz
	uint8_t pad[1];        // ensure 32 byte total size
};


struct BVHLittleTree {
	int startIndex;
	int numPrimitives;
	BVHNode* nodes;
	
};

struct BVH {
	BVH(std::vector<std::shared_ptr<Primitive>> p) : primitives(std::move(p)) {
		std::vector<BVHPrimitiveInfo> BVHPrimitives;
		BVHPrimitives.reserve(primitives.size());
		for (int i = 0; i < primitives.size(); i++) {
			BVHPrimitives.push_back({ i, primitives[i]->box });
		}
		MemoryArena arena(1024 * 1024);
		int totalNodes = 0;
		std::vector<std::shared_ptr<Primitive>> orderedPrimitives;
		orderedPrimitives.reserve(primitives.size());

		BVHNode* root;
		root = HLBVHBuild(arena, BVHPrimitives, &totalNodes, orderedPrimitives);
		primitives.swap(orderedPrimitives);
		BVHPrimitives.resize(0);
		printf("BVH created with %d nodes for %d "
			"primitives (%.4f MB), arena allocated %.2f MB\n",
			(int)totalNodes, (int)primitives.size(),
			float(totalNodes * sizeof(LinearBVHNode)) /
			(1024.f * 1024.f),
			float(arena.TotalAllocated()) /
			(1024.f * 1024.f));
		assert(root != nullptr);
		nodes = AllocAligned<LinearBVHNode>(totalNodes);
		int offset = 0;
		flattenBVHTree(root, &offset);

	}
	~BVH() { FreeAligned(nodes); }

	BVHNode* build(std::vector<MortonPrimitive>&, std::vector<Primitive>&);
	BVHNode* HLBVHBuild(MemoryArena& arena, const std::vector<BVHPrimitiveInfo>& BVHPrimitives, int* totalNodes, std::vector<std::shared_ptr<Primitive>>& orderedPrims);
	BVHNode* emit(BVHNode*& nodes, const std::vector<BVHPrimitiveInfo>& BVHPrimitives, MortonPrimitive* mortonPrimitives, std::vector<std::shared_ptr<Primitive>>&, int, int*, int*, int);
	BVHNode* buildSAH(MemoryArena& arena, std::vector<BVHNode*>& treeRoots, int start, int end, int* total) const;
	int flattenBVHTree(BVHNode*, int*);

	std::vector<std::shared_ptr<Primitive>> primitives;
	LinearBVHNode* nodes = nullptr;
	int maxPrimsInNode = 1;
};



inline uint32_t LeftShift3(uint32_t x) {
	if (x == (1 << 10)) --x;
	x = (x | (x << 16)) & 0b00000011000000000000000011111111;
	x = (x | (x << 8)) & 0b00000011000000001111000000001111;
	x = (x | (x << 4)) & 0b00000011000011000011000011000011;
	x = (x | (x << 2)) & 0b00001001001001001001001001001001;
	return x;
}

uint32_t EncodeMorton3(const Point3D& p) {
	return (LeftShift3(p.z) << 2) |
		   (LeftShift3(p.y) << 1) |
		   (LeftShift3(p.x) << 0);
}


short bitValue(uint32_t& number, uint32_t& mask) {
	return number & mask ? 1 : 0;
}


static void radixSort(std::vector<MortonPrimitive>* v)
{
	std::vector<MortonPrimitive> tempVector(v->size());
	const int bitsPerPass = 6;
	const int nBits = 30;
	static_assert((nBits % bitsPerPass) == 0,
		"Radix sort bitsPerPass must evenly divide nBits");
	const int nPasses = nBits / bitsPerPass;

	for (int pass = 0; pass < nPasses; ++pass) {
		// Perform one pass of radix sort, sorting _bitsPerPass_ bits
		int lowBit = pass * bitsPerPass;

		// Set in and out vector pointers for radix sort pass
		std::vector<MortonPrimitive>& in = (pass & 1) ? tempVector : *v;
		std::vector<MortonPrimitive>& out = (pass & 1) ? *v : tempVector;

		// Count number of zero bits in array for current radix sort bit
		const int nBuckets = 1 << bitsPerPass;
		int bucketCount[nBuckets] = { 0 };
		const int bitMask = (1 << bitsPerPass) - 1;
		for (const MortonPrimitive& mp : in) {
			int bucket = (mp.mortonCode >> lowBit) & bitMask;
			++bucketCount[bucket];
		}

		// Compute starting index in output array for each bucket
		int outIndex[nBuckets];
		outIndex[0] = 0;
		for (int i = 1; i < nBuckets; ++i)
			outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];

		// Store sorted values in output array
		for (const MortonPrimitive& mp : in) {
			int bucket = (mp.mortonCode >> lowBit) & bitMask;
			out[outIndex[bucket]++] = mp;
		}
	}
	// Copy final result from _tempVector_, if needed
	if (nPasses & 1) std::swap(*v, tempVector);
}

//BVHNode* BVH::build(std::vector<MortonPrimitive>& mortonPrimitives, std::vector<Primitive>& prims) {
//	
//	
//}

struct BucketInfo {
	int count = 0;
	BoundBox bounds;
};

BVHNode* BVH::HLBVHBuild(MemoryArena& arena, const std::vector<BVHPrimitiveInfo>& BVHPrimitives, int* totalNodes, std::vector<std::shared_ptr<Primitive>>& orderedPrims)  {
	BoundBox box;
	for (const BVHPrimitiveInfo& pi : BVHPrimitives) {
		box = box.Union(box, pi.centroid); // maybe it should be UNION @TODO
	}

	std::vector<MortonPrimitive> mortonPrims(BVHPrimitives.size());
	for (int i = 0; i < BVHPrimitives.size(); i++) {
		const int mortonBits = 10;
		const int mortonScale = 1 << mortonBits;

		mortonPrims[i].primitiveIndex = BVHPrimitives[i].primitiveNumber;
		Point3D p = box.offset(BVHPrimitives[i].centroid);
		p.x = p.x * mortonScale;
		p.y = p.y * mortonScale;
		p.z = p.z * mortonScale;
		mortonPrims[i].mortonCode = EncodeMorton3(p);
	}

	radixSort(&mortonPrims);

	//for (MortonPrimitive mp : mortonPrims) {
	//	std::cout << mp.primitiveIndex << " " << mp.mortonCode << std::endl;
	//}
	std::vector<BVHLittleTree> treesToBuild;

	uint32_t mask = 0b00111111111111000000000000000000; // first 12 bits describe the position of the primitive
	for (int start = 0, end = 1; end <= (int)mortonPrims.size(); end++) {
		if (end == mortonPrims.size() || ((mortonPrims[start].mortonCode & mask) != (mortonPrims[end].mortonCode & mask))) {
			int n = end - start;
			int maxNodes = 2 * n;
			BVHNode* nodes = arena.Alloc<BVHNode>(maxNodes, false);
			treesToBuild.push_back({ start, n, nodes });
			start = end;
		}
	}

	int orderedPrimsOffset = 0;
	orderedPrims.resize(primitives.size());
	int nodesCreated = 0;
	int firstBitIndex = 29 - 12;

	for (int i = 0; i < treesToBuild.size(); i++) {
		treesToBuild[i].nodes = BVH::emit(treesToBuild[i].nodes, BVHPrimitives, &mortonPrims[treesToBuild[i].startIndex], orderedPrims, treesToBuild[i].numPrimitives, &nodesCreated, &orderedPrimsOffset, firstBitIndex);
		*totalNodes += nodesCreated;
	}
	totalNodes += nodesCreated;
	std::vector<BVHNode*> finishedTrees;
	finishedTrees.reserve(treesToBuild.size());
	for (BVHLittleTree& tr : treesToBuild) {
		finishedTrees.emplace_back(tr.nodes);
	}
	return buildSAH(arena, finishedTrees, 0, finishedTrees.size(), totalNodes);

}

BVHNode* BVH::emit(BVHNode*& nodes, const std::vector<BVHPrimitiveInfo>& BVHPrimitive, MortonPrimitive* mortonPrimitives, std::vector<std::shared_ptr<Primitive>>& orderedPrimitives, int primitivesCount, int* totalNodes, int* orderedPrimsOffset, int bitIndex) {
	if (bitIndex == -1 || primitivesCount < maxPrimsInNode) {
		(*totalNodes)++;
		BVHNode* tmp = nodes++;
		BoundBox box;
		int firstPrimOffset = *orderedPrimsOffset;
		for (int i = 0; i < primitivesCount; i++) {
			int index = mortonPrimitives[i].primitiveIndex;
			orderedPrimitives[(int)(firstPrimOffset + i)] = primitives[index];
			box = box.Union(box, BVHPrimitive[index].box);
		}
		tmp->InitLeaf(0, primitivesCount, box);
		return tmp;

	}
	else {
		int mask = 1 << bitIndex;
		if ((mortonPrimitives[0].mortonCode & mask) == (mortonPrimitives[primitivesCount - 1].mortonCode & mask)){ // Next tree if nothing to split for this bit
			return emit(nodes, BVHPrimitive, mortonPrimitives, orderedPrimitives, primitivesCount, totalNodes, orderedPrimsOffset, bitIndex - 1);
		}
		int start = 0;
		int end = primitivesCount - 1;
		while (start + 1 != end) {
			int mid = (end - start) / 2 + start; // (start-end)/2
			if ((mortonPrimitives[start].mortonCode & mask) == (mortonPrimitives[mid].mortonCode & mask)) {
				start = mid;
			}
			else {
				end = mid;
			}
		}
		int split = end;
		(*totalNodes)++;
		BVHNode* tmp = nodes++;
		BVHNode* lbvh[2];
		lbvh[0] = emit(nodes, BVHPrimitive, mortonPrimitives, orderedPrimitives, split, totalNodes, orderedPrimsOffset, bitIndex-1);
		lbvh[1] = emit(nodes, BVHPrimitive, &mortonPrimitives[split], orderedPrimitives, primitivesCount - split, totalNodes, orderedPrimsOffset, bitIndex - 1);
		int axis = bitIndex % 3;
		tmp->InitInterior(axis, lbvh[0], lbvh[1]);
		return tmp;
	}
}

BVHNode* BVH::buildSAH(MemoryArena& arena, std::vector<BVHNode*>& treeRoots, int start, int end, int* total) const {
	
	int nodesCount = end - start;
	if (nodesCount == 1) {
		return treeRoots[start];
	}
	assert(nodesCount > 1);
	(*total)++;
	
	BVHNode* node = arena.Alloc<BVHNode>();
	BoundBox box;
	for (int i = start; i < end; i++) {
		box = Union(box, treeRoots[i]->box);
	}
	BoundBox centroidBox;
	for (int i = start; i < end; i++) {
		Point3D centroid = Point3D((treeRoots[i]->box.pMin.x + treeRoots[i]->box.pMax.x) * 0.5f, (treeRoots[i]->box.pMin.y + treeRoots[i]->box.pMax.y) * 0.5f, (treeRoots[i]->box.pMin.z + treeRoots[i]->box.pMax.z) * 0.5f);
		centroidBox = Union(centroidBox, centroid);
	}
	const int dimension = centroidBox.MaximumExtent() % 3;
	const int nBuckets = 12;
	struct Buckets {
		int count = 0;
		BoundBox box;
	};
	Buckets buckets[nBuckets];
	for (int i = start; i < end; i++) {
		float centroid = (treeRoots[i]->box.pMin[dimension] * 0.5f + treeRoots[i]->box.pMax[dimension] * 0.5f) ;
		int b = (int)(nBuckets * ((centroid - centroidBox.pMin[dimension]) / (centroidBox.pMax[dimension] - centroidBox.pMin[dimension])));
		if (b == nBuckets) b = nBuckets - 1;
		//assert(b < nBuckets);
		buckets[b].count++;
		buckets[b].box = Union(buckets[b].box, treeRoots[i]->box);
	}

	float cost[nBuckets - 1];
	for (int i = 0; i < nBuckets - 1; i++) {
		BoundBox b0, b1;
		int count0 = 0, count1 = 0;
		for (int j = 0; j <= i; j++) {
			b0 = Union(b0, buckets[j].box);
			count0 += buckets[j].count;	
		}
		for (int j = i+1; j < nBuckets; j++) {
			b1 = Union(b1, buckets[j].box);
			count1 += buckets[j].count;
		}

		cost[i] = (.125f + (count0 * b0.surfaceArea() + count1 * b1.surfaceArea())) / box.surfaceArea();
	}

	double minCost = cost[0];
	int minCostSplitBucket = 0;
	for (int i = 1; i < nBuckets - 1; ++i) {
		if (cost[i] < minCost) {
			minCost = cost[i];
			minCostSplitBucket = i;
		}
	}

	BVHNode** pmid = std::partition(&treeRoots[start], &treeRoots[(int)(end - 1)] + 1, [=](const BVHNode* node) {
			float centroid = (node->box.pMax[dimension] * 0.5f + node->box.pMin[dimension]*0.5f) ;
			int b = (int)(nBuckets * ((centroid - centroidBox.pMin[dimension]) / ((int)centroidBox.pMax[dimension] - centroidBox.pMin[dimension])));
			if (b == nBuckets) b = nBuckets - 1;
			return b <= minCostSplitBucket;
		});
	assert(pmid != nullptr);
	//std::cout << pmid << "  " << &treeRoots[0];
	int mid = (int)(pmid - &treeRoots[0]);
	//std::cout << start << " " << mid << std::endl;
	//std::cout << mid << " " << end << std::endl;
	std::cout << dimension << std::endl;
	assert(dimension < 3);

	node->InitInterior(dimension, this->buildSAH(arena, treeRoots, start, mid, total), this->buildSAH(arena, treeRoots, mid, end, total));
	return node;
}

int BVH::flattenBVHTree(BVHNode* node, int* offset) {
	LinearBVHNode* linearNode = &nodes[*offset];
	linearNode->bounds = node->box;
	int myOffset = (*offset)++;
	if (node->nPrimitives > 0) {
		linearNode->primitivesOffset = node->firstPrimOffset;
		linearNode->nPrimitives = node->nPrimitives;
	}
	else {
		// Create interior flattened BVH node
		linearNode->axis = node->splitAxis;
		linearNode->nPrimitives = 0;
		flattenBVHTree(node->children[0], offset);
		linearNode->secondChildOffset = flattenBVHTree(node->children[1], offset);
	}
	return myOffset;
}