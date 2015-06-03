FEWARPDIR = $(dir $(CURDIR))

# This needs to be set in order to find the febio .h files.
FEBINC = $(HOME)/febio-2.3.0/sdk/include
# This needs to be set in order to find the febio library files.
FEBLIB = $(HOME)/febio-2.3.0/sdk/lib

lnx64: PLAT = lnx64
lnx32: PLAT = lnx32
osx: PLAT = osx
gcc: PLAT = gcc
lnx64clean: PLAT = lnx64
lnx32clean: PLAT = lnx32
osxclean: PLAT = osx
gccclean: PLAT = gcc

export PLAT
export FEBLIB
export FEWARPDIR

# Assumes the Intel C++ compiler.  Use g++ for the gnu compiler and run "make gcc".
export CC = icpc
export FLG = -O3 -fPIC

# Can be removed if not using the Intel compiler.
export LNKFLG = -static-intel
export INC = -I$(FEBINC)

lnx64 lnx32 gcc:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="so" )

osx:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="dylib" )

lnx64clean lnx32clean gccclean:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="so" clean )

osxclean:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="dylib" clean )

.PHONY: lnx64 lnx32 osx gcc