add_executable(M2N main.c phasta.cc input.f readnblk.f cname.f setsyncioparam.cc new_interface.c dumbCvariables.c)
target_link_libraries(M2N phastaIO ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
        #force preprocessing
        set_source_files_properties(readnblk.f PROPERTIES COMPILE_FLAGS
"-qsuffix=cpp=f")
endif(CMAKE_Fortran_COMPILER_ID MATCHES XL)
if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
        #force preprocessing 
        set_source_files_properties(readnblk.f PROPERTIES COMPILE_FLAGS
"-Mpreprocess")
endif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
