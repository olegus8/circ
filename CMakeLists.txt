cmake_minimum_required(VERSION 3.1.0)
project(circ)

find_package(Boost 1.57.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  add_executable(circ circ.cpp)
endif()

