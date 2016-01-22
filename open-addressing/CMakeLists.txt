cmake_minimum_required(VERSION 2.8.11)
project(open-table)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 --coverage")

include_directories(../cpp-utils/include include)
set(SOURCE_FILES test/main.cpp include/open_table.hpp)
add_executable(open-table ${SOURCE_FILES})
