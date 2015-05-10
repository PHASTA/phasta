project(phastaIO)
enable_language(CXX)

# collect all the source files
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_Header(FCMangle.h)
include_directories("${phastaIO_BINARY_DIR}")
FILE(GLOB SRC_FILES *.cc *.f *.c)

ADD_LIBRARY(phastaIO ${SRC_FILES})

