# collect all the source files
FILE(GLOB SRC_FILES *.cc *.f *.c)

ADD_LIBRARY(phMetis ${SRC_FILES})

