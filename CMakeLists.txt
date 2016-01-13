cmake_minimum_required(VERSION 2.8.11)
project(Sandbox2)

message("System        : " ${CMAKE_SYSTEM})
message("CMake dir     : " ${CMAKE_ROOT})
message("Compiler      : " ${CMAKE_CXX_COMPILER})
message("Flags         : " ${CMAKE_CXX_FLAGS})

add_subdirectory(cpp-utils)
add_subdirectory(open-addressing)
add_subdirectory(perfect-hash)

