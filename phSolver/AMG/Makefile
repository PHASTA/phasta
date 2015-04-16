SUBSYSNAME := phasta/phSolver/AMG
MODULENAME := AMG
NVERS := 1.0
BUILDV := 1.0
NOSHARED = 1

PRECOMPILE := copymods

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
ifeq ($(ARCHOS),x86_linux-icc)
    DEFS := $(DEFS) -Dsun4_5
    FFLAGS := -w95 -W0 -cm -80
endif
ifeq ($(ARCHOS),alpha_osf)
    FFLAGS := -convert big_endian  -arch host -tune host -unroll 4 -fast -shared -w
    CCFLAGS := -convert big_endian -arch host -tune host -call_shared
    LDFLAGS := -ldxml -lsciport -lc -lfor -lFutil -lm -lots -lUfor
endif
ifeq ($(ARCHOS),ibm-bg)
    FFLAGS := $(FFLAGS) -qnullterm -qfixed=132
    DEFS := $(DEFS) -Dibm
endif

ifeq ($(ARCHOS),ibm_5)

    DEFS := $(DEFS) -Dibm6000
endif
ifeq ($(VERS),)
    DEFS := $(DEFS) -DLAUNCH_GDB
endif
ifeq ($(AMG),1)
    DEFS := $(DEFS) -DAMG
endif


dirs := .
DEPS := $(DEPS) 
mods  := ramg_data.f

INCLUDES := $(INCLUDES) -I../common  -I$(DEVROOT)/ARPACK/SRC
# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common

copymods:
	cp ../common/*.mod .

