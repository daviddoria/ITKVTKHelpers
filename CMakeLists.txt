cmake_minimum_required(VERSION 2.8)

project(ITKVTKHelpers)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

message("ITKVTKHelpers used from ${CMAKE_CURRENT_SOURCE_DIR}")

if(NOT TARGET ITKHelpers)
  add_subdirectory(ITKHelpers)
  include_directories(${ITKHelpers_includes})
  set(ITKVTKHelpers_includes ${ITKVTKHelpers_includes} ${ITKHelpers_includes})
  set(ITKVTKHelpers_libraries ${ITKVTKHelpers_libraries} ITKHelpers)
endif()

# Required for c++0x features used in Helpers sub-sub-module
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

# ITK
if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG ITKIOMeta
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology ITKDistanceMap)
  INCLUDE(${ITK_USE_FILE})
  set(ITKVTKHelpers_libraries ${ITKVTKHelpers_libraries} ${ITK_LIBRARIES})
endif()

# VTK
if(NOT VTK_FOUND)
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${VTK_USE_FILE})
  set(ITKVTKHelpers_libraries ${ITKVTKHelpers_libraries} ${VTK_LIBRARIES})
endif()

# Allow inheritance of include directories
if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(ITKVTKHelpers_includes ${ITKVTKHelpers_includes} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  set(ITKVTKHelpers_libraries ${ITKVTKHelpers_libraries} PARENT_SCOPE)
endif()

SET(ITKVTKHelpers_BuildTests OFF CACHE BOOL "ITKVTKHelpers_BuildTests")
if(ITKVTKHelpers_BuildTests)
  add_subdirectory(Tests)
endif(ITKVTKHelpers_BuildTests)

add_library(ITKVTKHelpers ITKVTKHelpers.cpp)
target_link_libraries(ITKVTKHelpers ${ITKVTKHelpers_libraries})
