cmake_minimum_required(VERSION 2.8)

project(ITKVTKHelpers)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# Only include ITKHelpers if it is not already available.
get_property(FoundITKHelpers GLOBAL PROPERTY ITKHelpersIncludeDirs SET)
if(NOT FoundITKHelpers)
  add_subdirectory(ITKHelpers)
endif()

get_property(ITKHelpersIncludeDirs GLOBAL PROPERTY ITKHelpersIncludeDirs)
set(ITKVTKHelpers_include_dirs ${ITKVTKHelpers_include_dirs} ${ITKHelpersIncludeDirs})
get_property(ITKHelpersLibs GLOBAL PROPERTY ITKHelpersLibs)
set(ITKVTKHelpers_libraries ${ITKVTKHelpers_libraries} ${ITKHelpersLibs})

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
endif()
set(ITKVTKHelpers_libraries ${ITKVTKHelpers_libraries} ${ITK_LIBRARIES})

# VTK
if(NOT VTK_FOUND)
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${VTK_USE_FILE})
endif()
set(ITKVTKHelpers_libraries ${ITKVTKHelpers_libraries} ${VTK_LIBRARIES})

add_custom_target(ITKVTKHelpersSources SOURCES
ITKVTKHelpers.h
ITKVTKHelpers.hpp
)

# Give the compiler all required include directories
include_directories(${ITKVTKHelpers_include_dirs})

# Create the library
add_library(ITKVTKHelpers ITKVTKHelpers.cpp)
target_link_libraries(ITKVTKHelpers ${ITKVTKHelpers_libraries})
set(ITKVTKHelpers_libraries ${ITKVTKHelpers_libraries} ITKVTKHelpers)

# Allow this project to be detected and used as a submodule
set(ITKVTKHelpers_include_dirs ${ITKVTKHelpers_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(GLOBAL PROPERTY ITKVTKHelpersIncludeDirs ${ITKVTKHelpers_include_dirs})
set_property(GLOBAL PROPERTY ITKVTKHelpersLibs ${ITKVTKHelpers_libraries})

# Build the tests if requested
option(ITKVTKHelpers_BuildTests "ITKVTKHelpers_BuildTests" OFF)
if(ITKVTKHelpers_BuildTests)
  add_subdirectory(Tests)
endif(ITKVTKHelpers_BuildTests)

# Display where this code was used from (if it is used as a submodule, there may be multiple instances of this submodule in the project, only the first of which is used)
option(ITKVTKHelpers_ShowSubmoduleLocation "Show the path from which ITKVTKHelpers was used?" OFF)
if(ITKVTKHelpers_ShowSubmoduleLocation)
  message("ITKVTKHelpers used from ${CMAKE_CURRENT_SOURCE_DIR}")
endif(ITKVTKHelpers_ShowSubmoduleLocation)
