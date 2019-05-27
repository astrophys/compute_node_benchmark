#!/bin/sh
# To Run : bash run_stream.sh /path/to/stream

#COREIDS=`egrep 'processor|physical id' /proc/cpuinfo  | \
#         perl -pe '/processor/ and s/\n/ /;' | \
#         egrep 'physical id.*: 0'| \
#         perl -ane 'print ("$F[2]\n");'`
set -e
COREIDS=`grep -P 'processor[\t ]' /proc/cpuinfo | cut -d: -f2 | tr -d ' '`
COUNT=1
STREAM=$1

for CORE in $COREIDS; do 
   if [ "$CORES" ] ; then
       CORES="$CORES,$CORE";
   else
       CORES=$CORE
   fi
   echo -n "$COUNT threads: "
   export OMP_NUM_THREADS=$COUNT
   taskset -c $CORES $STREAM | \
       egrep '^(Copy|Scale|Add|Triad):' | \
       awk '{print $1 " " $2}' | \
       tr '\n' ' '
    echo
    COUNT=`expr $COUNT + 1`
done
