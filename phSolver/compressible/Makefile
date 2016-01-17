SUBSYSNAME := phasta/phSolver/compressible
MODULENAME := compressible
BUILDV := 1.0
ISEXECUTABLE := 1
OUTNAME := phastaC.exe
ifeq ($(ARCHOS), )
    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

ifeq ($(ARCHOS),x86_cygwin)
  DEVROOT := $(shell cygpath -w -m $(DEVROOT))
  IOLIB_DIR := $(DEVROOT)/phasta/phastaIO/lib/$(ARCHOS)
  SHPLIB_DIR := $(DEVROOT)/phasta/shapeFunction/lib/$(ARCHOS)
  METIS_DIR := $(DEVROOT)/phasta/phMetis/lib/$(ARCHOS)
  WINIOLIB_DIR  := $(shell cygpath -w $(IOLIB_DIR) )
  WINSHPFN_DIR  := $(shell cygpath -w $(SHPLIB_DIR) )
  WINMETIS_DIR  := $(shell cygpath -w $(METIS_DIR) )
  INCLUDES := -I "C:\Program Files\MPICH\SDK\include"
  ADDSTDLIB := -LIBPATH:../lib/x86_cygwin libcommon.lib \
               -LIBPATH:"C:\Program Files\MPICH\SDK\Lib" mpich.lib \
               -LIBPATH:"$(WINIOLIB_DIR)" libphastaIO.lib \
               -LIBPATH:"$(WINSHPFN_DIR)" libshapeFunction.lib\
               -LIBPATH:"$(WINMETIS_DIR)" libmetis.lib\
                            Advapi32.lib Ws2_32.lib /subsystem:console /release 
else

  PRECOMPILE := copymods 
 ifeq ($(ARCHOS),x86_linux)
##    DEFS := $(DEFS) -Dsun4_5
    FFLAGS := -ffixed-form -ffixed-line-length-132
	LDFLAGS := $(LDFLAGS) -lgfortran -lgfortranbegin -lc -lm
 endif
 ifeq ($(ARCHOS),sun4_5)
       DEFS := $(DEFS) -Dsun4_5
       LDFLAGS := $(LDFLAGS) -lfsu -lfai2 -lfsumai -lfminlai -lfminvai \
                  -lfmaxvai -lfmaxlai -lF77 -lsunmath -lthread \
                  -lsocket -lnsl -lrt -lnsl -laio  -lsunperf -lmvec
  endif
  ifeq ($(ARCHOS),sgi_6-64)
       DEFS := $(DEFS) -Dsgi_6 -Dsgi 
       LDFLAGS := $(LDFLAGS) -lc -lfortran -lftn -lffio -lfpe -lffio -lfpe -lm
  endif
  ifeq ($(ARCHOS),ia64_linux)
        DEFS := $(DEFS) -Dsun4_5
  endif
  ifeq ($(ARCHOS),x86_linux-icc)
        DEFS := $(DEFS) -DLINUX 
        CFLAGS := $(CFLAGS) -w
        FFLAGS := $(FFLAGS) -w95 -W0 -cm -80
        LDFLAGS := $(LDFLAGS) -lifcore -ldl -lm -lpthread
  endif
  ifeq ($(ARCHOS),x86_64_linux-icc)
        DEFS := $(DEFS) -DLINUX 
        CFLAGS := $(CFLAGS) -w
        FFLAGS := $(FFLAGS) -w95 -W0 -cm -80
        LDFLAGS := $(LDFLAGS) -lifcore -ldl -lm -lpthread
  endif
  ifeq ($(ARCHOS),x86_64_linux-pgi)
        DEFS := $(DEFS) -DLINUX
        FFLAGS := $(FFLAGS) -Mextend
        LDFLAGS := $(LDFLAGS) -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl -lrt
#        LDFLAGS := $(LDFLAGS) -pgf90libs
  endif
  ifeq ($(ARCHOS),alpha_osf)
        DEFS := $(DEFS) -Ddecalp 
        FFLAGS := $(FFLAGS) -convert big_endian  -arch host -tune host -unroll 4\
                  -fast -shared -w
        CCFLAGS := $(CCFLAGS) -convert big_endian -arch host -tune host -call_shared
        LDFLAGS := $(LDFLAGS) -ldxml -lsciport -lc -lfor -lFutil -lm -lots -lUfor 
  endif
  ifeq ($(ARCHOS),ibm-bg)
        FFLAGS := $(FFLAGS) -qnullterm -qfixed=132
        DEFS := $(DEFS) -Dibm
        LDFLAGS := $(LDFLAGS) -L/opt/ibmcmp/xlf/bg/10.1/blrts_lib -lxlf90 -lxlfmath
  endif

endif

mods  := mod_timedata.f stats.f BCprofile.f mod_BCsfIDmap.f

dirs := .
DEPS := $(DEPS) phasta/phSolver/$(VERSION)/phSolver phasta/phastaIO phasta/phMetis\
        phasta/shapeFunction
LIBS := common phastaIO metis shapeFunction
ADDLIB = $(ADDSTDLIB)
INCLUDES := $(INCLUDES) -I../common

# need to include common for intel fortran .mod files
include $(DEVROOT)/Util/buildUtil/make.common

.PHONY: copymods
.PHONY: $(OUTFILE)

copymods:
	cp ../common/*.mod .
