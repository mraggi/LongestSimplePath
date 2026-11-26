cmake_minimum_required(VERSION 2.6)
project(lsp)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_definitions(-Wall)
add_executable(lsp main.cpp digraph.cpp pseudotopoorder.cpp helpers.cpp)
add_executable(lsp_unweighted main_read_from_stdin.cpp digraph.cpp pseudotopoorder.cpp helpers.cpp)
# add_executable(cumberbach main_cumberbach.cpp digraph.cpp pseudotopoorder.cpp helpers.cpp)
# add_executable(broughton main_broughton.cpp digraph.cpp pseudotopoorder.cpp helpers.cpp)
 
