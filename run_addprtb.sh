#!/bin/sh
#
# add rescaled perturbations to r_sig.f00
#
set -ex
#if [ $# -lt 2 ]; then
#  echo "Usage : ./run_addprtb.sh cycle bvstart"
#  exit 1
#fi
CYCLE=${1:-1}
BVSTART=${2:-1}
IDATE=${IDATE:-2022061800}
MEMBER=${MEMBER:-10}
BV_H=${INCCYCLE:-6}
TETYPE=${TETYPE:-dry}
SCL=${SCL}
QADJ=${QADJ:-yes} #super saturation and dry adjustment
ORTH=${ORTH:-yes} #orthogonalization
PSUB=${PSUB:-no}  #subtract perturbation
BP=${BP} #with boundary perturbation
head=${HEAD:-bv$TETYPE}
# directories
CDIR=` pwd `
SRCDIR=${CDIR}/build
DATA=${DATA:-${CDIR}/rsm2rsm27_bgmtest} #input top directory
OUT=${OUT:-${CDIR}/rsm2rsm27_bgm} #output top directory
BASE0=${BASE0:-${CDIR}/gfsp2rsm27_nomad} #(for the first cycle) arbitrary perturbation base directory
BASE1=${BASE1:-${CDIR}/gfsp2rsm27_rda} #(for the first cycle) arbitrary perturbation base directory

set +e
foffset=`expr ${BV_H} \* \( ${CYCLE} - 1 \)`
set -e
SDATE=`date -j -f "%Y%m%d%H" -v+${foffset}H +"%Y%m%d%H" "$IDATE"`
DATADIR=$DATA/$SDATE # input directory (must include unperturbed run)
OUTDIR=$OUT/$SDATE # output directory
if [ ! -d $DATADIR ]; then
  echo "No such directory : $DATADIR"
  exit 3
fi
EXEC=addprtb
cd $SRCDIR
make ${EXEC}
cd -
rm -rf tmp
mkdir -p tmp
cd tmp
rm -f *.grd
ln -fs ${SRCDIR}/${EXEC} ${EXEC}
# unperturbed run
echo $SDATE
ln -s $DATADIR/r_sig.f00 rb.0000.sig.grd
ln -s $DATADIR/r_sfc.f00 rb.0000.sfc.grd
# perturbation base
if [ do$PSUB = doyes ]; then
  SIGN=m
else
  SIGN=
fi
MEM=1
while [ $MEM -le $MEMBER ]; do
PMEM=`printf '%0.3d' $MEM` #prtb member
MEMDIR=$OUTDIR/${head}${SIGN}${PMEM}
## restart check
if [ -s $MEMDIR/r_sigi -a -s $MEMDIR/r_sigitdt -a -s $MEMDIR/r_sfci ]; then
  echo 'Restart files exist !!'
  MEM=`expr $MEM + 1`
#  continue
  exit
fi
mkdir -p $MEMDIR
##copy namelists
cp $DATADIR/rsmparm $MEMDIR/
cp $DATADIR/rsmlocation $MEMDIR/
cp $DATADIR/rfcstparm $MEMDIR/
cp $DATADIR/station.parm $MEMDIR/
##copy orography data
cp $DATADIR/rmtn.parm $MEMDIR/
cp $DATADIR/rmtnoss $MEMDIR/
cp $DATADIR/rmtnslm $MEMDIR/
cp $DATADIR/rmtnvar $MEMDIR/
MEM4=`printf '%0.4d' $MEM`
ln -s $MEMDIR/r_sig.f00 ro.$MEM4.sig.grd
ln -s $MEMDIR/r_sfc.f00 ro.$MEM4.sfc.grd
if [ $CYCLE -eq $BVSTART ]; then
## First cycle: interpolating perturbation base fields into the computational domain using rinp.sh
## In this demo script, we skip the interpolation procedure.
  MEM4=`expr 2 \* $MEM - 1`
  MEM4=`printf '%0.4d' $MEM4`
  cp $DATADIR/tmp_addprtb/ri.$MEM4.sig.grd .
#  cp $DATADIR/pdate.txt . # arbitrary chosen dates
#  irow=$MEM
#  PDATE1=`cat pdate.txt | awk '{if(NR == '$irow') {print $1}}'`
#  irow=`expr $irow + $MEMBER`
#  PDATE2=`cat pdate.txt | awk '{if(NR == '$irow') {print $1}}'`
#  echo $PDATE1 $PDATE2
#  BASE=$BASE0
#  PDATE=$PDATE1
#  if [ ! -d $BASE0/$PDATE ]; then
#    BASE=$BASE1
#  fi
#  if [ ! -d $BASE/$PDATE ]; then
#    echo "No such directory : $BASE/$PDATE"
#    exit 3
#  fi
#  mkdir -p $PDATE
#  cd $PDATE
#  hhr=12
#  cp $DATADIR/rsmlocation .
#  ln -s $DATADIR/rmtn.parm .
#  ln -s $DATADIR/rmtnoss .
#  ln -s $DATADIR/rmtnslm .
#  ln -s $DATADIR/rmtnvar .
#  if [ do$G2R = doyes ] ; then
#    ln -fs $BASE/$PDATE/sigf$hhr rb_sigf$hhr
#    ln -fs $BASE/$PDATE/sfcf$hhr rb_sfcf$hhr
#  fi
#  if [ do$P2R = doyes ] ; then
#    if [ do$CWBGFS = doyes ] ; then
#      if [ $hhr -lt 100 ] ; then hhrr=0$hhr ; fi
#      ln -fs $BASE/$PDATE/otgb2_$hhrr rb_pgbf$hhr
#    else
#      ln -fs $BASE/$PDATE/pgbf$hhr rb_pgbf$hhr
#    fi
#  else
#    if [ do$C2R = doyes ] ; then
#      cp $BASE/$PDATE/r_sig.f$hhr rb_sigf$hhr
#      cp $BASE/$PDATE/r_sfc.f$hhr rb_sfcf$hhr
#    fi
#  fi
#  #$USHDIR/rinp.sh $NEST $hhr > /dev/null 2>&1 || exit 10
#  $USHDIR/rinp.sh $NEST $hhr > rinp.log 2>&1 || exit 10
#  cp r_sigi ../ri.$MEM4.sig.grd
#  cd ..
  MEM4=`expr 2 \* $MEM`
  MEM4=`printf '%0.4d' $MEM4`
  cp $DATADIR/tmp_addprtb/ri.$MEM4.sig.grd .
#  #PDATE=`date -j -f "%Y%m%d%H" -v-24H +"%Y%m%d%H" "${PDATE}"`
#  BASE=$BASE0
#  PDATE=$PDATE2
#  if [ ! -d $BASE0/$PDATE ]; then
#    BASE=$BASE1
#  fi
#  if [ ! -d $BASE/$PDATE ]; then
#    echo "No such directory : $BASE/$PDATE"
#    exit 3
#  fi
#  mkdir -p $PDATE
#  cd $PDATE
#  #hhr=24
#  cp $DATADIR/rsmlocation .
#  ln -s $DATADIR/rmtn.parm .
#  ln -s $DATADIR/rmtnoss .
#  ln -s $DATADIR/rmtnslm .
#  ln -s $DATADIR/rmtnvar .
#  if [ do$G2R = doyes ] ; then
#    ln -fs $BASE/$PDATE/sigf$hhr rb_sigf$hhr
#    ln -fs $BASE/$PDATE/sfcf$hhr rb_sfcf$hhr
#  fi
#  if [ do$P2R = doyes ] ; then
#    if [ do$CWBGFS = doyes ] ; then
#      if [ $hhr -lt 100 ] ; then hhrr=0$hhr ; fi
#      ln -fs $BASE/$PDATE/otgb2_$hhrr rb_pgbf$hhr
#    else
#      ln -fs $BASE/$PDATE/pgbf$hhr rb_pgbf$hhr
#    fi
#  else
#    if [ do$C2R = doyes ] ; then
#      cp $BASE/$PDATE/r_sig.f$hhr rb_sigf$hhr
#      cp $BASE/$PDATE/r_sfc.f$hhr rb_sfcf$hhr
#    fi
#  fi
#  #$USHDIR/rinp.sh $NEST $hhr > /dev/null 2>&1 || exit 10
#  $USHDIR/rinp.sh $NEST $hhr > rinp.log 2>&1 || exit 10
#  cp r_sigi ../ri.$MEM4.sig.grd
#  cd ..
else
  # forecast difference
  fh2=$BV_H
  if [ $fh2 -lt 10 ]; then
    fh2=0$fh2
  fi
  PDATE=`date -j -f "%Y%m%d%H" -v-${BV_H}H +"%Y%m%d%H" "${SDATE}"` #a
  echo $PDATE #a
  MEM4=`expr 2 \* $MEM - 1`
  MEM4=`printf '%0.4d' $MEM4`
  ln -s $DATA/$PDATE/r_sig.f$fh2 ri.$MEM4.sig.grd
  MEM4=`expr 2 \* $MEM`
  MEM4=`printf '%0.4d' $MEM4`
  ln -s $DATA/$PDATE/${head}${PMEM}/r_sig.f$fh2 ri.$MEM4.sig.grd
fi
MEM=`expr $MEM + 1`
done #MEM -le MEMBER
### set namelist
if [ "$TETYPE" = "dry" ];then
  epsq=0.0d0
elif [ "$TETYPE" = "weak" ]; then
  epsq=0.1d0
else
  epsq=1.0d0
fi
if [ do$QADJ = doyes ];then
  adjust_q=T
fi
if [ do$ORTH = doyes ];then
  orth=T
fi
if [ do$PSUB = doyes ]; then
  lsub=T
fi
SPINUP=`expr 24 / $BV_H + 1`
if [ do$SCL != do ];then
  teref=${SCL}
  if echo "$SCL" | grep -q "^[0-9]\+$";then
    # SCL is integer
    teref=${SCL}.0d0
  fi
else
  teref=3.0d0
fi
#### set namelist
cat <<EOF >namelist
&namlst_prtb
 member=${MEMBER},
 setnorm=T,
 teref=${teref},
 epsq=${epsq},
! lonw=110.0,
! lone=153.0,
! lats=15.0,
! latn=47.0,
 adjust_q=${adjust_q},
 orth=${orth},
 lsub=${lsub},
&end
EOF
#
./${EXEC} < namelist 1>${EXEC}.log 2>${EXEC}.err
mv ${EXEC}.log ${EXEC}.err $DATADIR/
#
## prepare for the initial files of ensemble members
fh=00
MEM=1
while [ $MEM -le $MEMBER ]; do
PMEM=`printf '%0.3d' $MEM` #prtb member
MEMDIR=$OUTDIR/${head}${SIGN}${PMEM}
cp $MEMDIR/r_sig.f$fh $MEMDIR/r_sigi
cp $MEMDIR/r_sig.f$fh $MEMDIR/r_sigitdt
cp $MEMDIR/r_sfc.f$fh $MEMDIR/r_sfci
ls -l $MEMDIR
MEM=`expr $MEM + 1`
done
echo END
