#bring in the macros defined in the m4 directory.. these include AX_GMI, and AX_DEBUG
ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES        = libshapeFunction.la 
libshapeFunction_la_SOURCES  = \
 src/blend.c \
 src/blendDrv.c \
 src/Bn.c \
 src/BnDrv.c \
 src/En.c \
 src/EnDrv.c \
 src/Fn.c \
 src/FnDrv.c \
 src/Lagrange.c \
 src/modeShape.c \
 src/modeShapeDrv.c \
 src/parDrv.c \
 src/TriShapeAndDrv.cc \
 src/uniformP.c 

libshapeFunction_la_include_HEADERS = \
 src/shapeFuncInternals.h 

# The directory in which to install headers
libshapeFunction_la_includedir = $(includedir)

AM_CFLAGS= \
 -I$(top_srcdir) \
 $(DEBUG_FLAGS) \
 $(OPTIMIZATION_FLAGS) \
 $(PROFILING_FLAGS)

AM_CXXFLAGS=${AM_CFLAGS}
