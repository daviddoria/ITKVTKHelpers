# Allow headers in tests to be included like
# #include "ITKVTKHelpers.h" rather than needing
# #include "ITKVTKHelpers/ITKVTKHelpers.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(TestITKVTKHelpers TestITKVTKHelpers.cpp)
target_link_libraries(TestITKVTKHelpers ITKVTKHelpers)
