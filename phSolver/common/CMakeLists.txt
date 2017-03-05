# collect all the source files
FILE(GLOB SRC_FILES *.cc *.f *.f90 *.c *.F90)
set(PHASTA_COMMON_INCLUDE_DIRS)
set(PHASTA_COMMON_LIBS)
if(PHASTA_USE_MPI)
  add_definitions(-DOMPI_SKIP_MPICXX=1)
  set(PHASTA_COMMON_INCLUDE_DIRS ${PHASTA_COMMON_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${MPI_INCLUDE_PATH}
    ${MPI_Fortran_INCLUDE_PATH}
  )
  set(PHASTA_COMMON_LIBS ${PHASTA_COMMON_LIBS}
    ${MPI_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
  )
endif(PHASTA_USE_MPI)

if(NOT PHASTA_CHEF_ENABLED)
  set(SRC_FILES ${SRC_FILES} phstreamEmpty/phstreamEmpty.cc)
  set(PHASTA_COMMON_INCLUDE_DIRS
      ${PHASTA_COMMON_INCLUDE_DIRS}
      ${CMAKE_CURRENT_SOURCE_DIR}/phstreamEmpty)
endif()

include_directories(${PHASTA_COMMON_INCLUDE_DIRS})

ADD_LIBRARY(common ${SRC_FILES})
if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
        set_source_files_properties(fillsparse.f PROPERTIES COMPILE_FLAGS
"-qsuffix=cpp=f")
        set_source_files_properties(petscdbg.f PROPERTIES COMPILE_FLAGS
"-qsuffix=cpp=f")
endif(CMAKE_Fortran_COMPILER_ID MATCHES XL)
TARGET_LINK_LIBRARIES(common PUBLIC ${PHASTA_COMMON_LIBS})
if(PHASTA_CHEF_ENABLED)
  TARGET_LINK_LIBRARIES(common PUBLIC SCOREC::ph)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/phastaCommonConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/phastaCommonConfig.cmake"
)

get_git_head_revision(PHASTA_GIT_REFSPEC PHASTA_GIT_SHA1)
configure_file(phasta_version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/phasta_version.h )
set(git_hash_string "# git hash ${PHASTA_GIT_SHA1}")
execute_process(
  COMMAND echo ${git_hash_string}
  COMMAND cat - input.config
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_FILE ${CMAKE_BINARY_DIR}/input.config
)

if(PHASTA_TESTING)
  add_subdirectory(test)
endif(PHASTA_TESTING)
