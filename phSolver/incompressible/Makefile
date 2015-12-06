SUBSYSNAME := phasta/phSolver/incompressible
MODULENAME := incompressible
BUILDV := 1.0
ISEXECUTABLE := 1
OUTNAME := phastaIC.exe
ifeq ($(ARCHOS), )
    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

ifeq ($(ARCHOS),x86_cygwin)
# the /FORCE option is being used here because of IV8-MSVCRT.Lib 
# incompatabilities while including libles.lib (which was compiled on IV7).  
# It forces exec with missing symbols which could be dangerous but seems 
# to work so far

  LDFLAGS:= /FORCE /STACK:32000000
  DEVROOT := $(shell cygpath -w -m $(DEVROOT))
  DEFS = -DACUSIM_NT
  IOLIB_DIR := $(DEVROOT)/phasta/phastaIO/lib/$(ARCHOS)
  SHPLIB_DIR := $(DEVROOT)/phasta/shapeFunction/lib/$(ARCHOS)
  METIS_DIR := $(DEVROOT)/phasta/phMetis/lib/$(ARCHOS)
  WINLSLIB_DIR := $(shell cygpath -w $(LESLIBDIR)/lib/$(ARCHOS) )
  WINIOLIB_DIR  := $(shell cygpath -w $(IOLIB_DIR) )
  WINSHPFN_DIR  := $(shell cygpath -w $(SHPLIB_DIR) )
  WINMETIS_DIR  := $(shell cygpath -w $(METIS_DIR) )
  INCLUDES := -I "C:\Program Files\MPICH\SDK\include"
  ADDSTDLIB := -LIBPATH:../lib/x86_cygwin libcommon.lib \
               -LIBPATH:"C:\Program Files\MPICH\SDK\Lib" mpich.lib \
               -LIBPATH:"$(WINLSLIB_DIR)" libles.lib \
               -LIBPATH:"$(WINIOLIB_DIR)" libphastaIO.lib \
               -LIBPATH:"$(WINSHPFN_DIR)" libshapeFunction.lib\
               -LIBPATH:"$(WINMETIS_DIR)" libmetis.lib\
			    Advapi32.lib Ws2_32.lib /subsystem:console /release
else

  PRECOMPILE := copymods
  ARPACKPATH := $(DEVROOT)/ARPACK

  LIBARPACK := $(ARPACKPATH)/libarpack_linux.a
  LIBPARPACK := $(ARPACKPATH)/parpack_MPI-linux.a
  
  ifeq ($(ARCHOS),ia64_linux)
    DEFS := $(DEFS) -Dsun4_5 -DACUSIM_LINUXIPF
  endif
  ifeq ($(ARCHOS),sun4_5)
       DEFS := $(DEFS) -Dsun4_5 -DACUSIM_SUN
       LDFLAGS := $(LDFLAGS) -lfsu -lfai2 -lfsumai -lfminlai -lfminvai \
                  -lfmaxvai -lfmaxlai -lF77 -lsunmath -lthread \
                  -lsocket -lnsl -lrt -lnsl -laio  -lsunperf -lmvec
  endif
  ifeq ($(ARCHOS),sgi_6-64)
       DEFS := $(DEFS) -Dsgi_6 -Dsgi -DACUSIM_SGI64
       fortlib := /usr/lib64/libcomplib.sgimath.so
       LDFLAGS := $(LDFLAGS) -lc -lfortran $(fortlib) -lftn -lffio\
                  -lfpe -lffio -lfpe -lm
  endif 
  ifeq ($(ARCHOS),x86_linux-icc)
	  DEFS := $(DEFS) -DLINUX -DACUSIM_LINUX
	  CFLAGS := $(CFLAGS) -w
	  FFLAGS := $(FFLAGS) -w95 -W0 -cm -80 -fpp
	  LDFLAGS := $(LDFLAGS) -lifcore -lpthread
#          LDFLAGS := $(LDFLAGS) -L/usr/local/intel/mkl72/lib/32 -lmkl_ia32 -lmkl_p4 -lmkl_lapack
#	  LDFLAGS := (LDFLAGS) -lintrins -lCEPCF90 -lPEPCF90 -lF90 
  endif
  ifeq ($(ARCHOS),x86_64_linux-icc)
	  DEFS := $(DEFS) -DLINUX -DACUSIM_LINUX
	  CFLAGS := $(CFLAGS) -w
	  FFLAGS := $(FFLAGS) -w95 -W0 -cm -80 -fpp
	  LDFLAGS := $(LDFLAGS) -lifcore -lpthread
#	  LDFLAGS := $(LDFLAGS) -L/usr/local/intel/mkl72/lib/32 -lmkl_ia32 -lmkl_p4 -lmkl_lapack
#	  LDFLAGS := $(LDFLAGS) -lintrins -lCEPCF90 -lPEPCF90 -lF90
  endif
  ifeq ($(ARCHOS),x86_64_linux-pgi)
        DEFS := $(DEFS) -DLINUX
        FFLAGS := $(FFLAGS) -Mextend -Mpreprocess
        LDFLAGS := $(LDFLAGS) -lpgf90 -lpgf90_rpm1 -lpgf902 -lpgf90rtl -lpgftnrtl -lrt
#        LDFLAGS := $(LDFLAGS) -pgf90libs
  endif
  ifeq ($(ARCHOS),alpha_osf)
        DEFS := $(DEFS) -Ddecalp -DACUSIM_ALPHA
	FFLAGS := $(FFLAGS) -convert big_endian -arch host -tune host -unroll 4\
		  -fast -shared -w
	CCFLAGS := $(CCFLAGS) -convert big_endian -arch host -tune host -call_shared
	LDFLAGS := $(LDFLAGS) -ldxml -lsciport -lc -lfor -lFutil -lm -lots -lUfor -lmpi 
  endif
  ifeq ($(ARCHOS),ibm-bg)
        FFLAGS := $(FFLAGS) -qnullterm -qfixed=132 -qlist -qsource -qsuffix=cpp=f
        DEFS := $(DEFS) -Dibm
        LDFLAGS := $(LDFLAGS) -L/opt/ibmcmp/xlf/bg/10.1/blrts_lib -lxlf90 -lxlfmath
  endif
  ifeq ($(ARCHOS),ibm-bgp)
        FFLAGS := $(FFLAGS) -qnullterm -qfixed=132 -qlist -qsource -qsuffix=cpp=f
        DEFS := $(DEFS) -Dibm
        LDFLAGS := $(LDFLAGS) -L/soft/apps/ibmcmp/xlf/bg/11.1/bglib -lxlopt -lxl -lxlf90 -lxlfmath -L/opt/ibmcmp/xlsmp/bg/1.7/bglib -lxlomp_ser
  endif

endif

ifeq ($(AMG),1)
    DEFS := $(DEFS) -DAMG
    ifeq ($(ARCHOS),ibm-bg)
      FDEFS := $(FDEFS) '-WF,-DAMG'
    elif ($(ARCHOS),ibm-bgp)
      FDEFS := $(FDEFS) '-WF,-DAMG'
    else
      FDEFS := $(FDEFS) -DAMG
    endif
    LDFLAGS := $(LDFLAGS) $(LIBPARPACK) $(LIBARPACK)
endif


dirs := .
mods  := timedata.f stats.f
DEPS := $(DEPS) phasta/phSolver/$(VERSION)/phSolver phasta/phastaIO phasta/phMetis\
        phasta/shapeFunction LIBLES/1.5 
LIBS := $(LIBS) common phastaIO metis shapeFunction les

ifeq ($(AMG),1)
    DEPS := $(DEPS) phasta/phSolver/$(VERSION)/phSolver
    LIBS := $(LIBS) AMG
endif

INCLUDES := $(INCLUDES) -I../common
ADDLIB = $(ADDSTDLIB)
# need to include common for intel fortran .mod files
include $(DEVROOT)/Util/buildUtil/make.common

.PHONY: copymods
.PHONY: $(OUTFILE)

copymods:
	cp ../common/*.mod .
ifeq ($(AMG),1)	
	cp ../AMG/*.mod .
endif
