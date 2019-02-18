cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0048 NEW) # avoid annoing warnings
project(tvl1flow C)

message(STATUS "Configuring project: " ${PROJECT_NAME})

###############################################################################

#set(PROJECT_VERSION_MAJOR  0)
#set(PROJECT_VERSION_MINOR  1)
#set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
#    CACHE INTERNAL "The project version.")

###############################################################################

# NOTE: -fPIC is added so that the static library can be included into
#       shared libraries.
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -w -fPIC -O3 ")
set(CMAKE_BUILD_TYPE Release)

###############################################################################

# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else ()
	set (CMAKE_C_STANDARD 99)
endif ()

###############################################################################

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

###############################################################################

include_directories(${IIO_INCLUDE_DIR})

# tvl1 static library
add_library(tvl1flow STATIC tvl1flow_lib.c)
target_link_libraries(tvl1flow iio m)

# tvl1 executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
add_executable(tvl1flow-bin main.c)
target_link_libraries(tvl1flow-bin iio m)
set_target_properties(tvl1flow-bin PROPERTIES OUTPUT_NAME tvl1flow)

###############################################################################

# Expose public includes to other subprojects through cache variable
set(TVL1FLOW_INCLUDE_DIR ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)

