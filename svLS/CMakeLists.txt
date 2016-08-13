project(svLS Fortran)
cmake_minimum_required(VERSION 2.8)
if(NOT PHASTA_SKIP_MPI_CHECK)
find_package(MPI REQUIRED)
endif()
include_directories(${MPI_Fortran_INCLUDE_PATH})

SET(FSRCS ADDBCMUL.f
		BC.f
		#BICGS.f
		CGRAD.f
		COMMU.f
		CPUT.f
		DOT.f
		GE.f
		GMRES.f
		INCOMMU.f
		LHS.f
		LS.f
		NORM.f
		NSSOLVER.f
		#OMPLA.f
		PRECOND.f
		SOLVE.f
		SPARMUL.f)
add_library(svLS ${FSRCS})

install(TARGETS svLS
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

set(HEADERS
			svLS_API_c.h
			svLS_API.h
			svLS_c.h
			svLS.h
			svLS_STD.h
			svLS_STRUCT_c.h 
			svLS_STRUCT.h
	)

install(FILES ${HEADERS} DESTINATION include)

if(CMAKE_INSTALL_PREFIX)
set(libdir ${CMAKE_INSTALL_PREFIX}/lib)
set(includedir ${CMAKE_INSTALL_PREFIX}/include)
set(implicit_libs "")
foreach(L IN LISTS CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
	set(implicit_libs "${implicit_libs} ${CMAKE_LINK_LIBRARY_FLAG}${L}")
endforeach()
configure_file(svls.pc.in svls.pc)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/svls.pc DESTINATION lib/pkgconfig)
endif()
