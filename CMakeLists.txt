cmake_minimum_required(VERSION 2.8)

project(ITKVTKHelpers)

if(NOT TARGET ITKHelpers)
  add_subdirectory(ITKHelpers)
endif()

# Required for c++0x features used in Helpers sub-sub-module
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

SET(ITKVTKHelpers_BuildTests OFF CACHE BOOL "ITKVTKHelpers_BuildTests")
if(ITKVTKHelpers_BuildTests)
  add_subdirectory(Tests)
endif(ITKVTKHelpers_BuildTests)

add_library(ITKVTKHelpers ITKVTKHelpers.cpp)
target_link_libraries(ITKVTKHelpers ${ITK_LIBRARIES} ${VTK_LIBRARIES} ITKHelpers)
