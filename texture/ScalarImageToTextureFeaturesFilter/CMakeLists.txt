
cmake_minimum_required(VERSION 2.8)

project(ScalarImageToTextureFeaturesFilter)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()

add_executable(ScalarImageToTextureFeaturesFilter MACOSX_BUNDLE ScalarImageToTextureFeaturesFilter.cxx)

if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(ScalarImageToTextureFeaturesFilter ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(ScalarImageToTextureFeaturesFilter ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )

