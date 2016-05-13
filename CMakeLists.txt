cmake_minimum_required(VERSION 2.6)
project(lsp)

set(CMAKE_BUILD_TYPE Release)

SET (CMAKE_CXX_COMPILER "/usr/bin/clang++")
SET (CMAKE_AR      		"/usr/bin/llvm-ar")
SET (CMAKE_LINKER  		"/usr/bin/llvm-ld")
SET (CMAKE_NM      		"/usr/bin/llvm-nm")
SET (CMAKE_OBJDUMP 		"/usr/bin/llvm-objdump") 
SET (CMAKE_RANLIB  		"/usr/bin/llvm-ranlib")

add_definitions(-std=c++14)
add_definitions(-Ofast)

add_executable(lsp main.cpp moviegraph.cpp digraph.cpp path.cpp pseudotopoorder.cpp helpers.cpp trainer.cpp)
