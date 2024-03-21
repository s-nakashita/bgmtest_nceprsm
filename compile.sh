#!/bin/sh
## edit following your environment
CC=${CC:-gcc}
FC=${FC:-gfortran}
export CC FC

## If you have the precompiled w3lib, set the library locations.
W3LIBDIR=
W3LIB=
if [ -z $W3LIB ] || [ -z $W3LIBDIR ]; then
## compile w3lib
FFLAGS_4="-fno-range-check -fconvert=big-endian -fallow-argument-mismatch"
FFLAGS_8="-fdefault-integer-8 -fdefault-real-8 -fno-range-check -fconvert=big-endian -fallow-argument-mismatch"
FFLAGS_d="-fdefault-real-8 -fno-range-check -fconvert=big-endian -fallow-argument-mismatch"
ARFLAGS_R=" "
CFLAGS_R=" -g -DMACOS"
#CFLAGS_R=" -g -DLINUX" # comment out for Linux or FreeBSD
export FFLAGS_4 FFLAGS_8 FFLAGS_d ARFLAGS_R CFLAGS_R
cd lib
mkdir -p incmod/w3_4
mkdir -p incmod/w3_8
mkdir -p incmod/w3_d
cd src/w3lib-1.7
./makelibw3.sh
cd ../../..
fi

## compile addprtb & calcte
if [ -d build ]; then
rm -rf build
fi
if  [ -z $W3LIB ] || [ -z $W3LIBDIR ]; then
cmake -B build -DCMAKE_Fortran_COMPILER=${FC}
else
cmake -B build -DCMAKE_Fortran_COMPILER=${FC} -DW3LIBDIR=${W3LIBDIR} -DW3LIB=${W3LIB}
fi
cd build && make
