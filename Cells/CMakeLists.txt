cmake_minimum_required(VERSION 2.8.9)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_MACOSX_RPATH 1)

project (Cells_C)
set (SOURCE Cells_C.cpp)
add_library(Cells_C MODULE ${SOURCE})
