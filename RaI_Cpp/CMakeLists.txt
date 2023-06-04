cmake_minimum_required(VERSION 3.16.3)

project(Reliefs_as_Images)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})
set(ENABLE_PRECOMPILED_HEADERS OFF)

add_executable ( Reliefs_as_Images src/main.cpp )
target_link_libraries ( Reliefs_as_Images ${OpenCV_LIBS} ${CERES_LIBRARIES} )