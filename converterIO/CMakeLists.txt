project(converterIO)
cmake_minimum_required(VERSION 2.8.5)

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX=1)

include_directories(${MPI_INCLUDE_PATH})
add_executable(converterO2N converterO2N.cc)
add_executable(converterN2O converterN2O.cc)
target_link_libraries(converterO2N phastaIO)
target_link_libraries(converterN2O phastaIO)
