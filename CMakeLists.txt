cmake_minimum_required(VERSION 3.2)
project(RCM)

#set(CMAKE_OSX_ARCHITECTURES "x86_64")

option(YOCTO_OPENGL "Build OpenGL apps" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
if(NOT GENERATOR_IS_MULTI_CONFIG AND CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin/debug)
endif(NOT GENERATOR_IS_MULTI_CONFIG AND CMAKE_BUILD_TYPE MATCHES Debug)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
if(GENERATOR_IS_MULTI_CONFIG)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
endif(GENERATOR_IS_MULTI_CONFIG)

 include_directories(${CMAKE_SOURCE_DIR}/libs/eigen/)
#include_directories(/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3)

add_subdirectory(libs)
add_subdirectory(apps)


