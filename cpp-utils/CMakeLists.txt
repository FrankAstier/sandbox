cmake_minimum_required(VERSION 2.8) # Travis-ci doesn't want better
project(cpp-utils)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 --coverage")

#find_package(Boost 1.5.9 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(include)

set(SOURCE_FILES include/bytes.hpp include/fast_hash.hpp include/toolbox.hpp include/stl_io.hpp
        include/codecs.hpp include/type_traits.hpp)

add_executable(bytes-tests test/bytes_tests.cpp ${SOURCE_FILES})
add_executable(toolbox-tests test/toolbox_tests.cpp ${SOURCE_FILES})
add_executable(stl-io-tests test/stl_io_tests.cpp ${SOURCE_FILES})
add_executable(codecs-tests test/codecs_tests.cpp ${SOURCE_FILES})
add_executable(fast-hash-tests test/fast_hash_tests.cpp ${SOURCE_FILES})