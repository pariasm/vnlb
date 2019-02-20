cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0048 NEW) # avoid annoing warnings
project(VideoNLBayes C CXX Fortran)
message(STATUS "Configuring main project: " ${CMAKE_PROJECT_NAME})

###############################################################################

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ggdb")

set(CMAKE_BUILD_TYPE Release)

###############################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules/")
set(CMAKE_VERBOSE_MAKEFILE OFF)

###############################################################################

include_directories("$ENV{HOME}/local/include/")
set(CMAKE_PREFIX_PATH "$ENV{HOME}/local/")
message(STATUS "Setting CMAKE_PREFIX_PATH as " ${CMAKE_PREFIX_PATH} " for local libraries")
find_package( BLAS REQUIRED )
find_package( CBLAS REQUIRED )
find_package( LAPACKE REQUIRED )

###############################################################################

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

###############################################################################

add_subdirectory(lib/iio)
add_subdirectory(lib/tvl1flow)
add_subdirectory(src/VidUtils)
add_subdirectory(src/VNLBayes)
add_subdirectory(src)
add_subdirectory(scripts)
