#!/bin/bash
#SBATCH --partition=7302
#SBATCH --cpus-per-task=20

########### Get genomics data ##########

echo "Hello It Is RUN time!"
set -e
WORKING_PATH=`pwd`
echo ${WORKING_PATH}
NCORES=20
# export MKL_DEBUG_CPU_TYPE=5   # Only set for AMD cores

######## Matrix Multiply ########
# Set Environmental Variables
export PYTHONPATH=${HOME}/local/python/3.7/lib/python3.7/site-packages/
export PATH=${HOME}/local/python/3.7/bin:$PATH
export LD_LIBRARY_PATH=~/local/lib/:$LD_LIBRARY_PATH

export MKL_NUM_THREADS=1
export MKL_VERBOSE=1

for((CORE=0; CORE<${NCORES}; CORE++)); do
    python3 src/numpy_attempt.py & 
done
wait


echo "WARNING!!!" 
echo "  After inspecting the containts of ${WORKING_PATH}/benchmarking_out, you may "
echo "  consider deleteing it"
