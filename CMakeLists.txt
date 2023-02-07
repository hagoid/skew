cmake_minimum_required(VERSION 3.22)
project(skew)

set(CMAKE_CXX_STANDARD 17)

option(PROFILE_FLAG OFF)
if (PROFILE_FLAG)
    add_compile_definitions(PROFILE_FLAG)
endif(PROFILE_FLAG)

add_executable(skew skew.cpp)
add_executable(count_skew count_skew.cpp)
