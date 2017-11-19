cmake_minimum_required(VERSION 3.8)
project(MeanMax)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES player.cpp)
add_executable(MeanMax ${SOURCE_FILES})