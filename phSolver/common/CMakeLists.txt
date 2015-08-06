# collect all the source files
FILE(GLOB SRC_FILES *.cc *.f *.f90 *.c *.F90)
set(PHASTA_COMMON_INCLUDE_DIRS)
set(PHASTA_COMMON_LIBS)
if(PHASTA_USE_MPI)
  add_definitions(-DOMPI_SKIP_MPICXX=1)
  set(PHASTA_COMMON_INCLUDE_DIRS ${PHASTA_COMMON_INCLUDE_DIRS}
    ${MPI_INCLUDE_PATH}
    ${MPI_Fortran_INCLUDE_PATH}
  )
  set(PHASTA_COMMON_LIBS ${PHASTA_COMMON_LIBS}
    ${MPI_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
  )
endif(PHASTA_USE_MPI)

if(PHASTA_CHEF_ENABLED)
  set(PHASTA_COMMON_LIBS ${PHASTA_COMMON_LIBS} phstream)
  set(PHASTA_COMMON_INCLUDE_DIRS
    ${PHASTA_COMMON_INCLUDE_DIRS}
    ${PHASTACHEF_SOURCE_DIR}/phstream
  )
else()
  set(SRC_FILES ${SRC_FILES} phstreamEmpty/phstreamEmpty.cc)
  set(PHASTA_COMMON_INCLUDE_DIRS ${PHASTA_COMMON_INCLUDE_DIRS} phstreamEmpty)
endif()

include_directories(${PHASTA_COMMON_INCLUDE_DIRS})

ADD_LIBRARY(common ${SRC_FILES})
TARGET_LINK_LIBRARIES(common ${PHASTA_COMMON_LIBS})

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/phastaCommonConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/phastaCommonConfig.cmake"
)

if(PHASTA_TESTING)
  add_subdirectory(test)
endif(PHASTA_TESTING)
