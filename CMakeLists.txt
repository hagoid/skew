cmake_minimum_required(VERSION 3.22)
project(skew)

set(CMAKE_CXX_STANDARD 14)

add_executable(skew skew.cpp)
add_executable(count_skew count_skew.cpp)
