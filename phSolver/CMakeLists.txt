project(phSolver)
set(CMAKE_Fortran_MODULE_DIRECTORY ${phSolver_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)

set(CASES "/path/to/test/cases/"
  CACHE string
  "path to the test cases")
option(PHASTA_TESTING "enable tests" OFF)

option(PHASTA_USE_PETSC "Use PETSc GMRES Solver" OFF)
if(PHASTA_USE_PETSC)
find_package(PETSC REQUIRED)
link_directories(${PETSC_LIBRARY_DIRS})
include_directories(${PETSC_INCLUDE_DIRS})
add_definitions(-DHAVE_PETSC)
option(IS_JEDS_BRANCH "Use Jed's Experimental PETSc Optimizations" OFF)
if(IS_JEDS_BRANCH)
	add_definitions(-DJEDBROWN)
endif()
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

