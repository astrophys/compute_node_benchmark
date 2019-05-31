#!/bin/sh
## To run : bash run_kelvin.sh /path/to/outdir
OUTDIR=$1
INDIR=$2

if [ "$LD_LIBRARY_PATH" ] ; then
    export LD_LIBRARY_PATH=./:$LD_LIBRARY_PATH
else 
    export LD_LIBRARY_PATH=./
fi

CORES=`fgrep -c processor /proc/cpuinfo`
echo "CORES = $CORES"
echo "OUTDIR = $OUTDIR"
echo "INDIR  = $INDIR"
for N in `seq 1 $CORES` ; do 
    #kelvin ${INDIR}/kelvin.conf --PedigreeFile ${INDIR}/single.post > ${OUTDIR}/kelvin.out.$N 2>&1 &
    kelvin ${INDIR}/kelvin.conf > ${OUTDIR}/kelvin.out.$N 2>&1 &
done
wait

HOURS=0
MINS=0
SECS=0

for T in `fgrep 'Finished run' ${OUTDIR}/kelvin.out.* | cut -d' ' -f3 | tr -d '@,'`  ; do
    set -- `echo $T | tr 'hms' '   '`
    echo  "runtime string is $T, parses to $*"
    if [ $# -eq 3 ] ; then
        HOURS=`expr $HOURS + $1`
        MINS=`expr $MINS + $2`
        SECS=`expr $SECS + $3`
    elif [ $# -eq 2 ] ; then
        MINS=`expr $MINS + $1`
        SECS=`expr $SECS + $2`
    else
        SECS=`expr $SECS + $1`
    fi
done

echo hours $HOURS mins $MINS secs $SECS
HOURS=`expr $HOURS \* 3600`
SECS=`expr $SECS + $HOURS`
MINS=`expr $MINS \* 60`
SECS=`expr $SECS + $MINS`

echo Run time : `echo "scale = 4; print $SECS / $CORES" | bc`
