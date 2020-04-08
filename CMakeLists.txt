cmake_minimum_required(VERSION 3.15)
project(MortonCode_BVH_KD)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)



add_executable(MortonCode_BVH_KD
        geometry/BoundBox.hpp
        core/BVH.hpp
        core/memory.hpp
        geometry/Point3D.hpp
        core/Screen.hpp
    Source.cpp
        geometry/Vector3D.hpp
        )





