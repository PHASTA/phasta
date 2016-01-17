CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
#SET(CMAKE_FORTRAN_COMPILER gfortran gfortran)
PROJECT(PHASTA Fortran C CXX )
enable_language(Fortran)
enable_language(C)
enable_language(CXX)

enable_testing()
include(CTest)
set(MPIRUN "mpirun"
  CACHE string
  "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
  CACHE string
  "the command line flag to give process count to MPIRUN")

list(APPEND CMAKE_MODULE_PATH ${PHASTA_SOURCE_DIR}/CMakeFiles)

#Handle Fortran name mangling
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_Header(FCMangle.h SYMBOLS write_restart write_error write_displ write_field write_hessian)
include_directories("${PHASTA_BINARY_DIR}")

# Force static linking
IF(CMAKE_SYSTEM_NAME MATCHES BlueGeneQ-static)
  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
ENDIF(CMAKE_SYSTEM_NAME MATCHES BlueGeneQ-static)

#MESSAGE(ERROR ${CMAKE_Fortran_SOURCE_FILE_EXTENSIONS})
#SET(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f90 f F90 F)
# stuff for C/C++ calls to Fortran -- still needs to be done properly...
IF(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -ffixed-form -ffixed-line-length-132 -cpp " )
ENDIF(CMAKE_Fortran_COMPILER_ID MATCHES GNU)

IF(CMAKE_Fortran_COMPILER_ID MATCHES Intel)

SET(CMAKE_Fortran_FLAGS "-g -O3 -W0 -ip -80 -fpp ${CMAKE_Fortran_FLAGS}" )
message(STATUS "Found ifort")
message(STATUS ${CMAKE_Fortran_FLAGS})
ENDIF(CMAKE_Fortran_COMPILER_ID MATCHES Intel)

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mextend")
endif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
	#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qstrict -qnullterm -qfixed=132")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qnullterm -qfixed=132")
endif(CMAKE_Fortran_COMPILER_ID MATCHES XL)

IF(${CMAKE_BUILD_TYPE} MATCHES "Debug" )
  ADD_DEFINITIONS( -DDEBUG )
ENDIF(${CMAKE_BUILD_TYPE} MATCHES "Debug" )

# set archos "like" options for things like c calling fortran bindings
# later on we may want to do it dependent on behavior (e.g. 
# C_FORTRAN_CALL == uppercase, underline, nochange)
#OPTION(ARCHOS_OPTION "ArchOS dependent defines (e.g. -DLINUX)" )
SET(ARCHOS_OPTION "-DLINUX" CACHE STRING "ArchOS dependent defines (e.g. -DLINUX)" )
ADD_DEFINITIONS( ${ARCHOS_OPTION} )


SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  ${PHASTA_BINARY_DIR}/lib)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${PHASTA_BINARY_DIR}/bin)

OPTION(PHASTA_USE_MPI "Build MPI components" ON)
IF(PHASTA_USE_MPI)
  ADD_DEFINITIONS( -DPARALLEL -DMPI -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX )
  FIND_PACKAGE(MPI REQUIRED)
#  SET(CMAKE_MODULE_PATH ${PHASTA_SOURCE_DIR}/CMakeFiles)
  # bug in FIND_PACKAGE for Threads -- ignore for now (bug has been reported)...
  FIND_PACKAGE(Threads)
  IF(Threads_FOUND)
    INCLUDE(FindThreads)
  ENDIF(Threads_FOUND)
ENDIF(PHASTA_USE_MPI)

OPTION(PHASTA_BUILD_PHASTAIO "Build PhastaIO" ON)
if(PHASTA_BUILD_PHASTAIO)
	add_subdirectory(phastaIO)
	list(APPEND PHASTA_INCLUDES "${PHASTA_SOURCE_DIR}/phastaIO")
endif(PHASTA_BUILD_PHASTAIO)

OPTION(PHASTA_BUILD_SHAPEFUNCTION "Build shapeFunction" ON)
if(PHASTA_BUILD_SHAPEFUNCTION)
	list(APPEND PHASTA_INCLUDES "${PHASTA_SOURCE_DIR}/shapeFunction/src")
	add_subdirectory(shapeFunction)
endif(PHASTA_BUILD_SHAPEFUNCTION)

if(EXISTS ${PHASTA_SOURCE_DIR}/phShape)
OPTION(PHASTA_BUILD_PHSHAPE "Build phShape" OFF)
if(PHASTA_BUILD_PHSHAPE)
	list(APPEND PHASTA_INCLUDES "${PHASTA_SOURCE_DIR}/phShape/src")
	add_subdirectory(phShape)
endif(PHASTA_BUILD_PHSHAPE)
endif(EXISTS ${PHASTA_SOURCE_DIR}/phShape)

OPTION(PHASTA_USE_SVLS "Build svLS Linear Solver" OFF)
if(PHASTA_USE_SVLS)
	message(WARNING "\
svLS has a non-commercial license. \n \
by activating PHASTA_USE_SVLS in cmake \n \
you agree to the terms outlined in \
svLS/license.txt"
	)
	list(APPEND PHASTA_INCLUDES "${PHASTA_SOURCE_DIR}/svLS")
	add_subdirectory(svLS)
endif(PHASTA_USE_SVLS)

OPTION(PHASTA_BUILD_SOLVERS "Build Solvers" ON)
if(PHASTA_BUILD_SOLVERS)
	set(PHASTA_BUILD_PHASTAIO ON)
	if(PHASTA_USE_MPI)
		set(PHASTA_BUILD_METIS ON)
	endif(PHASTA_USE_MPI)
	list(APPEND PHASTA_INCLUDES "${PHASTA_SOURCE_DIR}/phSolver/common")
	include_directories(${PHASTA_INCLUDES})
	add_subdirectory(phSolver)
#TODO
endif(PHASTA_BUILD_SOLVERS)

OPTION(PHASTA_BUILD_CONVERTERIO "Build IO converter" ON)
if(PHASTA_BUILD_CONVERTERIO)
	add_subdirectory(converterIO)
endif(PHASTA_BUILD_CONVERTERIO)

OPTION(PHASTA_BUILD_ACUSTAT "Build AcuStat" ON)
if(PHASTA_BUILD_ACUSTAT)
	add_subdirectory(AcuStat)
endif(PHASTA_BUILD_ACUSTAT)

OPTION(PHASTA_BUILD_M2N "Build M2N" ON)
if(PHASTA_BUILD_M2N)
	add_subdirectory(M2N)
endif(PHASTA_BUILD_M2N)

OPTION(PHASTA_BUILD_M2NFixBnd "Build M2NFixBnd" ON)
if(PHASTA_BUILD_M2NFixBnd)
	add_subdirectory(M2NFixBnd)
endif(PHASTA_BUILD_M2NFixBnd)

if(EXISTS ${PHASTA_SOURCE_DIR}/phNSpre)
	add_subdirectory(phNSpre)
endif(EXISTS ${PHASTA_SOURCE_DIR}/phNSpre)

add_subdirectory(checkphasta)

