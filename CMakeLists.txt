cmake_minimum_required (VERSION 3.0)
project (ooccupy)
enable_language (CXX)
include_directories (./include)
file (GLOB SOURCES "src/*.cpp")
add_executable(ooccup++ ${SOURCES})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3  -std=c++17 -lstdc++fs -")
target_link_libraries(ooccup++  -lstdc++fs)
