project(svLS Fortran)
cmake_minimum_required(VERSION 2.8)
find_package(MPI REQUIRED)
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
