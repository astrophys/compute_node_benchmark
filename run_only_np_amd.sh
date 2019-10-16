#!/bin/bash
#SBATCH --partition=7302
#SBATCH --cpus-per-task=32

########### Get genomics data ##########

echo "Hello It Is RUN time!"
set -e
WORKING_PATH=`pwd`
echo ${WORKING_PATH}
NCORES=32
export MKL_DEBUG_CPU_TYPE=5

######## Matrix Multiply ########
# Set Environmental Variables
export PYTHONPATH=${HOME}/local/python/3.7/lib/python3.7/site-packages/
export PATH=${HOME}/local/python/3.7/bin:$PATH
export LD_LIBRARY_PATH=~/local/lib/:$LD_LIBRARY_PATH
rm -rf ${WORKING_PATH}/benchmarking_out_numpy/

for((CORE=0; CORE<${NCORES}; CORE++)); do
    mkdir -p ${WORKING_PATH}/benchmarking_out_numpy/${CORE}/data/matrix
    python3 src/numpy_attempt.py build_mat_mult_data ${WORKING_PATH}/benchmarking_out_numpy/${CORE}/ & 
done
wait


echo "WARNING!!!" 
echo "  After inspecting the containts of ${WORKING_PATH}/benchmarking_out, you may "
echo "  consider deleteing it"
