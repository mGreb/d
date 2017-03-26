cmake_minimum_required(VERSION 3.6)
project(cs_days)

FIND_PACKAGE(ALGLIB REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(cs_days ${SOURCE_FILES})