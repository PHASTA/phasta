add_executable(checkphasta checkphasta.cpp)
target_link_libraries(checkphasta common phastaIO ${MPI_LIBRARIES})
