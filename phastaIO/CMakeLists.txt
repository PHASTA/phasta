project(phastaIO)
enable_language(CXX)

find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_PATH})
# collect all the source files
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_Header(FCMangle.h)
include_directories("${phastaIO_BINARY_DIR}")
add_definitions(-DMPICH_SKIP_MPICXX)
add_definitions(-DOMPI_SKIP_MPICXX=1)
set(SRC_FILES phiotmrc.cc phastaIO.cc)

ADD_LIBRARY(phastaIO ${SRC_FILES})
target_link_libraries(phastaIO ${MPI_C_LIBRARIES})
