add_executable(checkphasta checkphasta.cpp)
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(checkphasta common phastaIO ${MPI_LIBRARIES})
