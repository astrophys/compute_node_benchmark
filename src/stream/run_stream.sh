#!/bin/sh
# To Run : bash run_stream.sh /path/to/stream

#COREIDS=`egrep 'processor|physical id' /proc/cpuinfo  | \
#         perl -pe '/processor/ and s/\n/ /;' | \
#         egrep 'physical id.*: 0'| \
#         perl -ane 'print ("$F[2]\n");'`
set -e
STREAM=$1
THREADS=1
NODES=

SOCKETIDS=`lscpu -e | grep -v SOCKET | tr -s ' ' | cut -d' ' -f3 | sort -u`
for SOCKET in $SOCKETIDS ; do
    NODEIDS=`lscpu -e | grep -v SOCKET | tr -s ' ' | cut -d' ' -f2,3 | egrep $SOCKET'$' | cut -d' ' -f1 | sort -u`
    for NODE in $NODEIDS ; do
	if [ "$NODES" ] ; then
	    NODES="$NODES,$NODE"
        else
            NODES=$NODE
        fi
	CPUIDS=`lscpu -e | grep -v SOCKET | tr -s ' ' | cut -d' ' -f1,2,3 | egrep "$NODE $SOCKET\$" | cut -d' ' -f1 | sort -nu`
	for CPU in $CPUIDS ; do
            echo -n THREADS $THREADS NODES $NODES CPU $CPU" "
	    OMP_NUM_THREADS=$THREADS numactl --membind=$NODES --cpunodebind=$NODES ${STREAM} | \
                egrep '^(Copy|Scale|Add|Triad):' | \
                awk '{print $1 " " $2}' | \
                tr '\n' ' '
            echo
	    THREADS=`expr $THREADS + 1`
        done
	
    done
done


exit 0
