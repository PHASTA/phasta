add_executable(checkphasta checkphasta.cpp)
target_link_libraries(checkphasta phastaIO ${MPI_LIBRARIES})
