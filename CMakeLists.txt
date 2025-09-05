cmake_minimum_required(VERSION 3.15)
project(GFN2_XTB_Project CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/gfn2_xtb" CACHE PATH "Install path" FORCE)

# Global optimization
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -fopenmp")

# Subdirectories
add_subdirectory(lib)
add_subdirectory(programs)

