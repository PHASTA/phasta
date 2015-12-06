project(phSolver)
set(CMAKE_Fortran_MODULE_DIRECTORY ${phSolver_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

set(CASES "/path/to/test/cases/"
  CACHE string
  "path to the test cases")
option(PHASTA_TESTING "enable tests" OFF)

option(PHASTA_USE_PETSC "Use PETSc GMRES Solver" OFF)
if(PHASTA_USE_PETSC)
find_package(PETSc REQUIRED)
set (BUILD_SHARED_LIBS NO) #FIXME - the petsc cmake config
                           #FIXME - file forces shared libs on...
                           #FIXME - which causes linking problems for us
find_path(PETSC_INC petscsys.h HINTS /usr/include ${PETSC_PACKAGE_INCLUDES})
find_path(PETSC_LIB libpetsc.a HINTS /usr/lib ${PETSC_PACKAGE_INCLUDES}../lib)
include_directories(${PETSC_PACKAGE_INCLUDES} ${PETSC_INC})
link_directories(${PETSC_LIB})
add_definitions(-DHAVE_PETSC)
endif(PHASTA_USE_PETSC)

add_subdirectory(common)
option(PHASTA_INCOMPRESSIBLE "Build the incompressible solver" OFF)
option(PHASTA_COMPRESSIBLE "Build the compressible solver" ON)
if(PHASTA_INCOMPRESSIBLE)
add_subdirectory(incompressible)
endif(PHASTA_INCOMPRESSIBLE)
if(PHASTA_COMPRESSIBLE)
add_subdirectory(compressible)
endif(PHASTA_COMPRESSIBLE)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/phastaConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/phastaConfig.cmake"
)

