cmake_minimum_required(VERSION 2.8)
project(SplineEditor)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#####################
# libigl
#####################

# libigl options: choose between header only and compiled static library
# Header-only is preferred for small projects. For larger projects the static build
# considerably reduces the compilation times
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" OFF)

# libigl options: choose your dependencies (by default everything is OFF except opengl) 
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"           ON)
option(LIBIGL_WITH_CGAL        "Use CGAL"           ON)
if(LIBIGL_WITH_CGAL) # Do not remove or move this block, the cgal build system fails without it
  find_package(CGAL REQUIRED)
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
  include(${CGAL_USE_FILE})
endif()

find_package(LIBIGL QUIET)
if (NOT LIBIGL_FOUND)
   message(FATAL_ERROR "libigl not found --- You can download it using: \n git clone --recursive https://github.com/libigl/libigl.git ${PROJECT_SOURCE_DIR}/../libigl")
endif()

set(LIBIGL_WITH_NANOGUI ON)
set(LIBIGL_WITH_OPENGL  ON)
set(LIBIGL_WITH_GLFW    ON)
set(LIBIGL_WITH_VIEWER  ON)

# Compilation flags: adapt to your needs 
if(MSVC)
  # Enable parallel compilation
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") 
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  # Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

# Adding libigl: choose the path to your local copy libigl 
# This is going to compile everything you requested 
add_subdirectory("${LIBIGL_INCLUDE_DIR}/../shared/cmake" "libigl")

# libigl information 
message("libigl includes: ${LIBIGL_INCLUDE_DIRS}")
message("libigl libraries: ${LIBIGL_LIBRARIES}")
message("libigl extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("libigl extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("libigl definitions: ${LIBIGL_DEFINITIONS}")

include_directories(${LIBIGL_INCLUDE_DIRS})
add_definitions(${LIBIGL_DEFINITIONS})

#####################
include_directories(src)
SET(TINYSPLINE_SRC src/tinyspline/tinysplinecpp.cpp src/tinyspline/tinyspline.c)
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fpermissive" )
#####################

FILE(GLOB SRCFILES src/*.cpp)

add_executable(${PROJECT_NAME} ${SRCFILES} ${TINYSPLINE_SRC} ${LIBIGL_EXTRA_SOURCES})
target_link_libraries(${PROJECT_NAME} ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES} alglib)

