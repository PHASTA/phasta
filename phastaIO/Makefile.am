#bring in the macros defined in the m4 directory.. these include AX_GMI, and AX_DEBUG
ACLOCAL_AMFLAGS = -I m4

SUBDIRS= .

lib_LTLIBRARIES        = libphastaio.la 
libphastaio_la_SOURCES  = \
phastaIO.cc

libphastaio_la_include_HEADERS = \
phastaIO.h

# The directory in which to install headers
libphastaio_la_includedir = $(includedir)

AM_CXXFLAGS= \
        -I$(top_srcdir) \
	$(DEBUG_FLAGS) \
    	$(OPTIMIZATION_FLAGS) \
    	$(PROFILING_FLAGS)
