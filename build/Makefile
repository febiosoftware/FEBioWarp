FEWARPDIR = $(dir $(CURDIR))

# This needs to be set in order to find the febio .h files.
FEBINC = $(HOME)/FEBio2
# This needs to be set in order to find the febio library files.
FEBLIB = $(HOME)/FEBio2/build/lib

lnx64d: PLAT = lnx64d
lnx64: PLAT = lnx64
lnx32: PLAT = lnx32
osx: PLAT = osx
osxd: PLAT = osxd
gcc: PLAT = gcc
lnx64clean: PLAT = lnx64
lnx64dclean: PLAT = lnx64d
lnx32clean: PLAT = lnx32
osxclean: PLAT = osx
osxdclean: PLAT = osxd
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

lnx64 lnx64d lnx32 gcc:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="so" )

osx osxd:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="dylib" )

lnx64clean lnx64dclean lnx32clean gccclean:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="so" clean )

osxclean osxdclean:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="dylib" clean )

.PHONY: lnx64 lnx64d lnx32 osx osxd gcc