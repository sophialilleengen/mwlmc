cmake_minimum_required(VERSION 3.14)
project(
  MWLMC
  VERSION 2.0.0
  HOMEPAGE_URL https://github.com/sophialilleengen/mwlmc
  LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using intel
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qno-offload")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()


### the next line might be necessary for Intel machines or MX machines 
### that built Python with the Interl infrastucture
### but will lead to problems with native Mac MX Python installations
### uncomment if neccessary
##set(CMAKE_OSX_ARCHITECTURES "x86_64")

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# required compiler features
#find_package(pybind11)

#find_package(Python COMPONENTS Interpreter Development)

find_package(Eigen3 REQUIRED)
add_subdirectory(extern/pybind11-2.9.2)

if(pybind11_FOUND)
  set(HAVEPYBIND11 TRUE)
  message(STATUS "Do we have pybind11? ${HAVEPYBIND11}")
endif()

if(Eigen3_FOUND)
  set(HAVEEIGEN TRUE)
  message(STATUS "Do we have eigen? ${HAVEEIGEN}")
  include_directories( ${EIGEN3_INCLUDE_DIRS} )
endif()

# add_subdirectory(src)
# add_subdirectory(examples)

set(common_LINKLIB ${EIGEN3_LIBRARIES})

# is this the best way to do this?
include_directories(include)

if(HAVEPYBIND11)
  pybind11_add_module(model src/mwlmc.cc )
endif()
