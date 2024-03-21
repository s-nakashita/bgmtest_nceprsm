#!/bin/sh
###############################################################
#
#   AUTHOR:    Gilbert - W/NP11
#
#   DATE:      01/11/1999
#
#   PURPOSE:   This script uses the make utility to update the libw3 
#              archive libraries.
#              It first reads a list of source files in the library and
#              then generates a makefile used to update the archive
#              libraries.  The make command is then executed for each
#              archive library, where the archive library name and 
#              compilation flags are passed to the makefile through 
#              environment variables.
#
#   REMARKS:   Only source files that have been modified since the last
#              library update are recompiled and replaced in the object
#              archive libraries.  The make utility determines this
#              from the file modification times.
#
#              New source files are also compiled and added to the object 
#              archive libraries.
#
###############################################################

#
#     Generate a list of object files that corresponds to the
#     list of Fortran ( .f ) files in the current directory
#
for i in `ls *.f`
do
  obj=`basename $i .f`
  OBJS="$OBJS ${obj}.o"
done
#
#     Generate a list of object files that corresponds to the
#     list of C ( .c ) files in the current directory
#
for i in `ls *.c`
do
  obj=`basename $i .c`
  OBJS="$OBJS ${obj}.o"
done
#
#     Remove make file, if it exists.  May need a new make file
#     with an updated object file list.
#
if [ -f make.libw3 ] 
then
  rm -f make.libw3
fi
#
#     Generate a new make file ( make.libw3), with the updated object list,
#     from this HERE file.
#
cat > make.libw3 << EOF
SHELL=/bin/sh

\$(LIB):	\$(LIB)( ${OBJS} )

.f.a:
	xlf -c \$(FFLAGS) \$<
	ar -ruv \$(AFLAGS) \$@ \$*.o
	#rm -f \$*.o

.c.a:
	/usr/vac/bin/xlc -c \$(CFLAGS) \$<
	ar -ruv  \$(AFLAGS) \$@ \$*.o
	#rm -f \$*.o
EOF
#
#     Update 4-byte version of libw3_4.a
#
export LIB="../../libw3_4.a"
export FFLAGS=" -O3 -qnosave -qmoddir=../../incmod/w3_4 -I../../incmod/w3_4"
export AFLAGS=" -X64"
export CFLAGS=" -O3 -q64 -DIBM4"
make -f make.libw3
#
#     Update 8-byte version of libw3_8.a
#
export LIB="../../libw3_8.a"
export FFLAGS=" -O3 -qnosave -qintsize=8 -qrealsize=8 -qmoddir=../../incmod/w3_8 -I../../incmod/w3_4"
export AFLAGS=" -X64"
export CFLAGS=" -O3 -q64 -DIBM4"
make -f make.libw3
#
#     Update Double Precision (Size of Real 8-byte and default Integer) version 
#     of libw3_d.a
#
export LIB="../../libw3_d.a"
export FFLAGS=" -O3 -qnosave -qrealsize=8 -qmoddir=../../incmod/w3_d -I../../incmod/w3_4"
export AFLAGS=" -X64"
export CFLAGS=" -O3 -q64 -DIBM4"
make -f make.libw3

rm -f make.libw3
