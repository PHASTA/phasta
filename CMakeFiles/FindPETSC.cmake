find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
	pkg_check_modules(PETSCPKG "PETSc")
endif(PKG_CONFIG_FOUND)

if(PETSCPKG_FOUND AND (NOT PETSc_DIR)) 
#if PETSc_DIR we're probably doing things the old way
#so just skip to that
#otherwise, try and use pkg-config
set(PETSC_LIBRARY_DIRS ${PETSCPKG_LIBRARY_DIRS})
set(PETSC_INCLUDE_DIRS ${PETSCPKG_INCLUDE_DIRS})
set(PETSC_LIBRARIES ${PETSCPKG_LIBRARIES})
set(PETSC_FOUND TRUE)

else()
find_package(PETSc REQUIRED)
if(PETSc_FOUND)
	set(PETSC_FOUND TRUE)
endif(PETSc_FOUND)

find_path(PETSC_INC petscsys.h HINTS /usr/include ${PETSC_PACKAGE_INCLUDES})
find_path(PETSC_LIB libpetsc.a HINTS /usr/lib ${PETSC_PACKAGE_INCLUDES}../lib)
set(PETSC_INCLUDE_DIRS ${PETSC_INC})
set(PETSC_LIBRARY_DIRS ${PETSC_LIB})
set(PETSC_LIBRARIES ${PETSC_LIB}/libpetsc.a ${PETSC_PACKAGE_LIBS})
endif(PETSCPKG_FOUND AND (NOT PETSc_DIR))
