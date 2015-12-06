SUBSYSNAME := phasta/phSolver/common
MODULENAME := common
NVERS := 1.0
BUILDV := 1.0
NOSHARED = 1
ifeq ($(ARCHOS), )
    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch )
endif

ifeq ($(ARCHOS),x86_cygwin)
    DEVROOT := $(shell cygpath -w -m $(DEVROOT))
    INCLUDES := $(INCLUDES) -I "C:\Program Files\MPICH\SDK\include"
    DEFS := $(DEFS) -Dintel
endif
ifeq ($(ARCHOS),sun4_5)
    DEFS := $(DEFS) -Dsun4_5 
endif
ifeq ($(ARCHOS),sgi_6-64)
    DEFS := $(DEFS) -Dsgi_6 -Dsgi
endif
ifeq ($(ARCHOS),ia64_linux)
    DEFS := $(DEFS) -Dsun4_5
endif
ifeq ($(ARCHOS),x86_linux)
    DEFS := $(DEFS) -DLINUX
    FFLAGS := -ffixed-form -ffixed-line-length-132
endif
ifeq ($(ARCHOS),x86_linux-icc)
    DEFS := $(DEFS) -DLINUX
    FFLAGS := $(FFLAGS) -w95 -W0 -cm -80
endif
ifeq ($(ARCHOS),x86_linux-icc-amd32)
    DEFS := $(DEFS) -DLINUX
    FFLAGS := $(FFLAGS) -w95 -W0 -cm -80
endif
ifeq ($(ARCHOS),x86_64_linux-icc)
    DEFS := $(DEFS) -DLINUX
    FFLAGS := $(FFLAGS) -w95 -W0 -cm -80
endif
ifeq ($(ARCHOS),x86_64_linux-pgi)
    DEFS := $(DEFS) -DLINUX
    FFLAGS := $(FFLAGS) -Mextend
endif
ifeq ($(ARCHOS),alpha_osf)
    FFLAGS := $(FFLAGS) -convert big_endian  -arch host -tune host -unroll 4 -fast -shared -w
    CCFLAGS := $(CCFLAGS) -convert big_endian -arch host -tune host -call_shared
    LDFLAGS := $(LDFLAGS) -ldxml -lsciport -lc -lfor -lFutil -lm -lots -lUfor
endif
ifeq ($(ARCHOS),ibm-bg)
    FFLAGS := $(FFLAGS) -qnullterm -qfixed=132
    DEFS := $(DEFS) -Dibm
endif
ifeq ($(ARCHOS),ibm-bgp)
    FFLAGS := $(FFLAGS) -qnullterm -qfixed=132
    DEFS := $(DEFS) -Dibm
endif
ifeq ($(VERS),)
    DEFS := $(DEFS) -DLAUNCH_GDB
endif

ifeq ($(AMG),1)
    DEFS := $(DEFS) -DAMG
endif

dirs := .
DEPS := $(DEPS) phasta/phastaIO 
mods  := spebc.f pointer.f \
        turbsa.f \
        perprep.f local_mass.f \
        filtprep.f  pvsqbi.f \
        dtn.f bardmc.f \
        turbke.f\
        bctint.f readnblk.f\
        mod_slpw.f

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common
