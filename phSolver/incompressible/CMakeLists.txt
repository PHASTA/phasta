file(GLOB SRC_FILES *.cc *.f *.f90 *.c)
list(REMOVE_ITEM SRC_FILES "main.cc")
add_library(incompressible ${SRC_FILES})
find_library(ACUSOLVE_LIB libles)
if(PHASTA_USE_MPI)
	include_directories(${MPI_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
endif(PHASTA_USE_MPI)
include_directories(${PHASTA_BINARY_DIR})
add_executable(phastaIC.exe main.cc)
set_target_properties(phastaIC.exe PROPERTIES HAS_CXX TRUE)
set_target_properties(incompressible PROPERTIES HAS_CXX TRUE)
set_target_properties(phastaIC.exe PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(incompressible PROPERTIES LINKER_LANGUAGE Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
	#force preprocessing of solfar.f
	set_source_files_properties(solfar.f PROPERTIES COMPILE_FLAGS "-qsuffix=cpp=f")
endif(CMAKE_Fortran_COMPILER_ID MATCHES XL)
if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
	#force preprocessing of solfar.f
	set_source_files_properties(solfar.f PROPERTIES COMPILE_FLAGS "-Mpreprocess")
endif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
#incompressible reuses a fortran module from common w/o telling cmake
#this is why we specify a module directory in the next cmakelists up
#and we need to serialize incompressible after common
add_dependencies(incompressible common)
if(PHASTA_USE_MPI)
	target_link_libraries(phastaIC.exe ${MPI_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}) 
endif(PHASTA_USE_MPI)

target_link_libraries(phastaIC.exe ${ACUSOLVE_LIB})
target_link_libraries(phastaIC.exe common incompressible common phastaIO)
if(PHASTA_BUILD_SHAPEFUNCTION)
	target_link_libraries(phastaIC.exe shapeFunction)
endif(PHASTA_BUILD_SHAPEFUNCTION)
if(PHASTA_BUILD_PHSHAPE)
	target_link_libraries(phastaIC.exe phshape)
endif(PHASTA_BUILD_PHSHAPE)
target_link_libraries(phastaIC.exe ${ACUSOLVE_LIB} incompressible)

if(PHASTA_TESTING)
  add_subdirectory(test)
endif(PHASTA_TESTING)
