cmake_minimum_required(VERSION 2.8.11)
project(linear_hash)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -g") # --coverage")

include_directories(../cpp-utils/include include)

set(SOURCE_FILES test/main.cpp
        include/linear-hash-LH.hpp
        include/linear-hash-LH1.hpp
        include/linear-hash-PS1.hpp
        include/linear-hash-PS2.hpp)

add_executable(linear-hash ${SOURCE_FILES})
