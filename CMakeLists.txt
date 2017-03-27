cmake_minimum_required(VERSION 2.6)
project(lsp)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_definitions(-Wall)
add_executable(lsp main.cpp digraph.cpp pseudotopoorder.cpp helpers.cpp)
 
