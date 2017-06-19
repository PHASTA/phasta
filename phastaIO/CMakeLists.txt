project(phastaIO)
enable_language(CXX)


if(NOT PHASTA_SKIP_MPI_CHECK)
find_package(MPI REQUIRED)
endif()
include_directories(${MPI_C_INCLUDE_PATH})
# collect all the source files
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_Header(FCMangle.h)
include_directories("${phastaIO_BINARY_DIR}")
add_definitions(-DMPICH_SKIP_MPICXX)
add_definitions(-DOMPI_SKIP_MPICXX=1)
set(SRC_FILES phiotmrc.cc phastaIO.cc phiompi.c)
if(NOT PHASTA_CHEF_ENABLED)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/phSolver/modules)
  set(SRC_FILES ${SRC_FILES}
    phiotimer_empty/phiotimer_empty.cc
    phiotimer_empty/phiotimer_empty.f)
  include_directories(phiotimer_empty)
endif()

ADD_LIBRARY(phastaIO ${SRC_FILES})
target_link_libraries(phastaIO PRIVATE ${MPI_C_LIBRARIES})
if(PHASTA_CHEF_ENABLED)
  TARGET_LINK_LIBRARIES(phastaIO PUBLIC SCOREC::ph)
endif()
