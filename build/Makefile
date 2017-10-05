FEWARPDIR = $(dir $(CURDIR))

# This needs to be set in order to find the febio .h files.
FEBINC = $(HOME)/FEBio2
# This needs to be set in order to find the febio library files.
FEBLIB = $(HOME)/FEBio2/build/lib

lnx64d: PLAT = lnx64d
lnx64s: PLAT = lnx64s
lnx64: PLAT = lnx64
lnx32: PLAT = lnx32
osx: PLAT = osx
osxd: PLAT = osxd
osxs: PLAT = osxs
gcc: PLAT = gcc
lnx64clean: PLAT = lnx64
lnx64dclean: PLAT = lnx64d
lnx64sclean: PLAT = lnx64s
lnx32clean: PLAT = lnx32
osxclean: PLAT = osx
osxdclean: PLAT = osxd
osxsclean: PLAT = osxs
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

lnx64 lnx64d lnx64s lnx32 gcc:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="so" )

osx osxd osxs:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="dylib" )

lnx64clean lnx64dclean lnx64sclean lnx32clean gccclean:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="so" clean )

osxclean osxdclean osxsclean:
	( cd $(PLAT); $(MAKE) -f ../fewarp.mk SFX="dylib" clean )

.PHONY: lnx64 lnx64d lnx64s lnx32 osx osxd osxs gcc