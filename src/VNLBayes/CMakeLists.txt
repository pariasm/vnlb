set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
add_library(vnlb SHARED VideoNLBayes.cpp LibMatrix.cpp)
target_link_libraries(vnlb vidutils)
install(TARGETS vnlb LIBRARY DESTINATION lib COMPONENT libs)
