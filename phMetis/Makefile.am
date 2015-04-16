#bring in the macros defined in the m4 directory.. these include AX_GMI, and AX_DEBUG
ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES        = libmetis.la 
libmetis_la_SOURCES  = \
 Lib/balance.c \
 Lib/bucketsort.c \
 Lib/ccgraph.c \
 Lib/coarsen.c \
 Lib/compress.c \
 Lib/debug.c \
 Lib/estmem.c \
 Lib/fm.c \
 Lib/fortran.c \
 Lib/frename.c \
 Lib/graph.c \
 Lib/initpart.c \
 Lib/kmetis.c \
 Lib/kvmetis.c \
 Lib/kwayfm.c \
 Lib/kwayrefine.c \
 Lib/kwayvolfm.c \
 Lib/kwayvolrefine.c \
 Lib/match.c \
 Lib/mbalance2.c \
 Lib/mbalance.c \
 Lib/mcoarsen.c \
 Lib/memory.c \
 Lib/mesh.c \
 Lib/meshpart.c \
 Lib/mfm2.c \
 Lib/mfm.c \
 Lib/mincover.c \
 Lib/minitpart2.c \
 Lib/minitpart.c \
 Lib/mkmetis.c \
 Lib/mkwayfmh.c \
 Lib/mkwayrefine.c \
 Lib/mmatch.c \
 Lib/mmd.c \
 Lib/mpmetis.c \
 Lib/mrefine2.c \
 Lib/mrefine.c \
 Lib/mutil.c \
 Lib/myqsort.c \
 Lib/ometis.c \
 Lib/parmetis.c \
 Lib/pmetis.c \
 Lib/pqueue.c \
 Lib/refine.c \
 Lib/separator.c \
 Lib/sfm.c \
 Lib/srefine.c \
 Lib/stat.c \
 Lib/subdomains.c \
 Lib/timing.c \
 Lib/util.c 

libmetis_la_include_HEADERS = \
 Lib/defs.h \
 Lib/macros.h \
 Lib/metis.h \
 Lib/proto.h \
 Lib/struct.h \
 Lib/rename.h 

# The directory in which to install headers
libmetis_la_includedir = $(includedir)

AM_CFLAGS= \
 -I$(top_srcdir) \
 $(DEBUG_FLAGS) \
 $(OPTIMIZATION_FLAGS) \
 $(PROFILING_FLAGS)
