#!/bin/sh
#
# calculate perturbation moist total energy from unperturbed run
#
set -ex
export OMP_NUM_THREADS=1
#if [ $# -lt 2 ]; then
#  echo "Usage : ./run_calcte.sh cycle"
#  exit 1
#fi
CYCLE=${1:-1}
MEMBER=${MEMBER:-10}
IDATE=${IDATE:-2022061800}
BV_H=${BV_H:-6}
TETYPE=${TETYPE:-dry}
SCL=${SCL}
QADJ=${QADJ:-no} #super saturation and dry adjustment
ORTH=${ORTH:-no}
BP=${BP} #with boundary perturbation
SCLBASE=${SCLBASE}
SCLPOW=${SCLPOW}

# directories
CDIR=` pwd `
SRCDIR=${CDIR}/build
DATADIR=${CDIR}/rsm2rsm27_bgmtest
if [ ! -d $DATADIR ]; then
  echo "No such directory : $DATADIR"
  exit 3
fi
WORKDIR=tmp
EXEC=calcte
cd $SRCDIR
make ${EXEC}
cd -
mkdir -p $WORKDIR
cd $WORKDIR
rm -f fort.*
ln -fs ${SRCDIR}/${EXEC} ${EXEC}
if [ $CYCLE -gt 1 ];then
  fhs=`expr $BV_H \* \( $CYCLE - 1 \)`
else
  fhs=0
fi
CDATE=`date -j -f "%Y%m%d%H" -v+${fhs}H +"%Y%m%d%H" "${IDATE}"` #a
dte=$BV_H
inch=1

MEM=1
while [ $MEM -le $MEMBER ];do
MEM=`printf '%0.3d' $MEM`
if [ $CYCLE -gt 1 ] && [ $BV_H -gt 6 ];then
  WDIR=bv${TETYPE}${SCL}${BV_H}h${MEM}${BP}${SCLBASE}
else
  WDIR=bv${TETYPE}${SCL}${MEM}${BP}${SCLBASE}
fi
if [ ! -z $SCLPOW ]; then
  WDIR=${WDIR}p${SCLPOW}
fi
if [ do$QADJ = doyes ];then
  WDIR=${WDIR}_qadj
fi
if [ do$ORTH = doyes ];then
  WDIR=${WDIR}_orth
fi

for dt in $(seq 0 $inch $dte);do
  fh=$dt
  if [ $fh -lt 10 ]; then
    fh=0$fh
  fi
# unperturbed run
  nsig=11
  ln -s $DATADIR/$CDATE/r_sig.f$fh fort.$nsig
# ensemble member
  nsig=`expr $nsig + 1`
  ln -s $DATADIR/$CDATE/${WDIR}/r_sig.f$fh fort.$nsig
cat <<EOF >NAMELIST
&NAMLST_PRTB
 lprtb=T,
 epsq=,
 kmax=,
 lonw=,
 lone=,
 lats=,
 latn=,
&END
EOF
  ./${EXEC} < NAMELIST #1>>${EXEC}.log 2>&1
  cat te.dat
  mv te.dat $DATADIR/$CDATE/${WDIR}/te${fh}h.dat #c
  rm fort.*
done #dt
MEM=`expr $MEM + 1`
done #MEMBER
echo END
