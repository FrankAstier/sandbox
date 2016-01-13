cmake_minimum_required(VERSION 2.8.11)
project(perfect_hash)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11") # --coverage")

message("System        : " ${CMAKE_SYSTEM})
message("CMake dir     : " ${CMAKE_ROOT})
message("Compiler      : " ${CMAKE_CXX_COMPILER})
message("Flags         : " ${CMAKE_CXX_FLAGS})

include_directories(cpp-utils/include)
set(SOURCE_FILES test/main.cpp)
add_executable(perfect-hash ${SOURCE_FILES})
