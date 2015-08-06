file(GLOB SRC_FILES *.cc *.f *.f90 *.c)
list(REMOVE_ITEM SRC_FILES "main.cc")
find_package(phastaCommon REQUIRED PATHS ${CMAKE_BINARY_DIR})
set(PHASTAC_LIBS ${CMAKE_THREAD_LIBS_INIT} ${PHASTA_COMMON_LIBS})

add_library(compressible ${SRC_FILES})

add_executable(phastaC.exe main.cc)
include_directories(
  ${PHASTAC_INCLUDE_DIRS}
  ${PHASTA_COMMON_INCLUDE_DIRS}
  ${PHASTA_BINARY_DIR}
)
set_target_properties(phastaC.exe PROPERTIES HAS_CXX TRUE)
set_target_properties(compressible PROPERTIES HAS_CXX TRUE)
set_target_properties(phastaC.exe PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(compressible PROPERTIES LINKER_LANGUAGE Fortran)
set(PHASTAC_LIBS ${PHASTAC_LIBS} phastaIO)
if(PHASTA_BUILD_SHAPEFUNCTION)
  set(PHASTAC_LIBS ${PHASTAC_LIBS} shapeFunction)
endif(PHASTA_BUILD_SHAPEFUNCTION)
if(PHASTA_BUILD_PHSHAPE)
  set(PHASTAC_LIBS ${PHASTAC_LIBS} phshape)
endif(PHASTA_BUILD_PHSHAPE)
#ughhhhh
target_link_libraries(phastaC.exe
  compressible ${PHASTAC_LIBS}
  compressible ${PHASTAC_LIBS}
)

if(PHASTA_TESTING)
  add_subdirectory(test)
endif(PHASTA_TESTING)
