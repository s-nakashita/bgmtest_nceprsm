#!/bin/sh
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
export OBJS
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

\$(LIB):	\$(OBJS)
	ar -ruv \$(ARFLAGS) \$@ *.o

.c.a:
	\$(CC) -c \$(CFLAGS) \$<

.f.a:
	\$(FC) -c \$(FFLAGS) \$<
EOF
#
#     Update 4-byte version of libw3_4.a
#
export LIB="../../libw3_4.a"
export FFLAGS=$FFLAGS_4
export CFLAGS=$CFLAGS_R
export ARFLAGS=$ARFLAGS_R
make -f make.libw3
mv *.mod ../../incmod/w3_4/.
rm -f *.o
#
#     Update 8-byte version of libw3_8.a
#
export LIB="../../libw3_8.a"
export FFLAGS=$FFLAGS_8
export CFLAGS=$CFLAGS_R
export ARFLAGS=$ARFLAGS_R
make -f make.libw3
mv *.mod ../../incmod/w3_8/.
rm -f *.o
#
#     Update Double Precision (Size of Real 8-byte and default Integer) version 
#     of libw3_d.a
#
export LIB="../../libw3_d.a"
export FFLAGS=$FFLAGS_d
export CFLAGS=$CFLAGS_R
export ARFLAGS=$ARFLAGS_R
make -f make.libw3
mv *.mod ../../incmod/w3_d/.
rm -f *.o
#
rm -f make.libw3
