set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_executable(vnlbayes main_vnlb.cpp)
add_executable(awgn awgn.cpp)
add_executable(psnr psnr.cpp)

target_link_libraries(vnlbayes vnlb ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(awgn vidutils)
target_link_libraries(psnr vidutils)

# These are not necessary
# ${CBLAS_LIBRARIES}
# ${IIO_LIBRARIES})

## external excludes
#include_directories(${LAPACK_INCLUDES})
#include_directories(${CBLAS_INCLUDE_DIR})
#include_directories(${IIO_INCLUDE_DIR})

