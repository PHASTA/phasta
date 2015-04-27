
# collect all the source files
FILE(GLOB SRC_FILES *.cc *.f *.f90 *.c *.F90)
#FIND_PACKAGE(MPI REQUIRED)   
if(PHASTA_USE_MPI)
	add_definitions(-DOMPI_SKIP_MPICXX=1)
	include_directories(${MPI_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})
endif(PHASTA_USE_MPI)
ADD_LIBRARY(common ${SRC_FILES})
if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
        set_source_files_properties(fillsparse.f PROPERTIES COMPILE_FLAGS
"-qsuffix=cpp=f")
        set_source_files_properties(petscdbg.f PROPERTIES COMPILE_FLAGS
"-qsuffix=cpp=f")
endif(CMAKE_Fortran_COMPILER_ID MATCHES XL)
TARGET_LINK_LIBRARIES(common ${MPI_LIBRARIES})
TARGET_LINK_LIBRARIES(common ${MPI_Fortran_LIBRARIES})
