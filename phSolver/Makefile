#
# in case the user wants to give the executables distinct names , he/she can specify the variable
# UXT on the make command line.  For instance 
# gmake VERS=opt COMPRESSIBLE=1 UXT=-myrinet will create a target
# phasta-compressible-myrinet-O
# 
#******NOTE: user has to preappend a '-' to the actual extension.
#

BUILDDIRS=common

ifeq ($(AMG),1)
  BUILDDIRS += AMG
endif

ifeq ($(COMPRESSIBLE),1)
  BUILDDIRS += compressible
else 
  BUILDDIRS += incompressible
endif

all:
	for i in ${BUILDDIRS}; do ( \
	  cd $$i; \
	  $(MAKE)) ; done

setup:
	for i in ${BUILDDIRS}; do ( \
	  cd $$i; \
	  $(MAKE) setup ) ; done

clean:
	for i in ${BUILDDIRS}; do ( \
	  cd $$i; \
	  $(MAKE) clean ) ; done
