cmake_minimum_required(VERSION 2.8.11)
project(linear_hash)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")

include_directories(../cpp-utils/include)
set(SOURCE_FILES test/main.cpp)
add_executable(linear-hash ${SOURCE_FILES})
