cmake_minimum_required(VERSION 3.5)

SET(CMAKE_VERBOSE_MAKEFILE "ON")
SET(CMAKE_COLOR_MAKEFILE "ON")

project(cs_days)

FIND_PACKAGE(ALGLIB REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -Wall -Wextra -Wunused")

set(SOURCE_FILES src/main.cpp src/P3D_thomas.cpp)
add_executable(cs_days ${SOURCE_FILES})
target_link_libraries(cs_days alglib)
