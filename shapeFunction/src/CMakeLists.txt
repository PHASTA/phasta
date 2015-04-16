# collect all the source files
FILE(GLOB SRC_FILES *.cc *.f *.c)

ADD_LIBRARY(shapeFunction ${SRC_FILES})
