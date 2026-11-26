cmake_minimum_required(VERSION 3.20)
project(longestsimplepath LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(LSP_BUILD_TESTS "Build unit tests (needs GoogleTest)" ON)

add_library(lsp
  src/ap.cpp src/graph.cpp src/io.cpp src/structure.cpp src/ordering.cpp src/dfs.cpp src/pto.cpp src/solver.cpp)
target_include_directories(lsp PUBLIC include)
target_compile_options(lsp PRIVATE -Wall -Wextra -Wpedantic)

add_executable(lsp_cli apps/lsp.cpp)
set_target_properties(lsp_cli PROPERTIES OUTPUT_NAME lsp)
target_link_libraries(lsp_cli PRIVATE lsp)

add_executable(lsp_check apps/lsp_check.cpp)
target_link_libraries(lsp_check PRIVATE lsp)

add_executable(lsp_example examples/minimal.cpp)
target_link_libraries(lsp_example PRIVATE lsp)

if(LSP_BUILD_TESTS)
  find_package(GTest)
  if(GTest_FOUND)
    enable_testing()
    add_executable(lsp_tests tests/test_main.cpp)
    target_link_libraries(lsp_tests PRIVATE lsp GTest::gtest GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(lsp_tests)
  else()
    message(WARNING "GoogleTest not found: tests disabled")
  endif()
endif()
