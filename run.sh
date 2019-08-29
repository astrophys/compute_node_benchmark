#!/bin/bash
# Author : Ali Snedden
# Date : 8/29/19
# Purpose:
#   This will run all the code that the Singularity container does. 
#   However, it will require that you've previously run the build.sh script and 
#   that you've got the extraneous/dependent software installed yourself
#
#   This script expects to be run within the within the compute_node_benchmark/
#   directory. You'll need to ensure that PYTHONPATH and PATH have been modified
#   appropriately to find pytho3 and Numpy. 
#
# HOW TO RUN (with SLURM) : 
#   sbatch --cpus-per-task=XX --mail-type=FAIL,REQUEUE,TIME_LIMIT_90 --mail-user=XX@YY.com run.sh
#
#
#

echo "Hello It Is RUN time!"
set -e
WORKING_PATH=`pwd`

######## Matrix Multiply ########
# Set Environmental Variables
#export PYTHONPATH=${HOME}/local/python/3.7/lib/python3.7/site-packages/
#export PATH=${HOME}/local/python/3.7/bin:$PATH
export PATH=${WORKING_PATH}/software/tophat-2.1.1:$PATH
export PATH=${WORKING_PATH}/software/hisat2-2.1.0:$PATH
export PATH=${WORKING_PATH}/software/bowtie2-2.3.5.1:$PATH
export PATH=${WORKING_PATH}/software/cufflinks-2.2.1:$PATH
export PATH=${WORKING_PATH}/software/samtools-1.9/bin:$PATH
export PATH=${WORKING_PATH}/software/stringtie-1.3.6:$PATH

if [ -d "benchmarking_out" ]; then 
    rm -r ${WORKING_PATH}/benchmarking_out
fi
mkdir -p ${WORKING_PATH}/benchmarking_out
mkdir -p ${WORKING_PATH}/benchmarking_out/data
mkdir -p ${WORKING_PATH}/benchmarking_out/output

python3 src/driver.py build_mat_mult_data ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py mat_mult_cache_opt ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py mat_mult_non_cache_opt ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py build_rnaseq_data ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py align_rnaseq_tophat ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py align_rnaseq_hisat ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py cufflinks_assemble ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py cuffmerge ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py cuffcompare ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py cuffquant ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py cuffnorm ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py cuffdiff ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/
python3 src/driver.py kelvin ${WORKING_PATH}/benchmarking_out/ ${WORKING_PATH}/ref/


### Run stream multiple times ###
echo "Running STREAM_ARRAY_SIZE = 10M" 
echo "Run 1"
bash src/stream/run_stream.sh src/stream/stream.10M
echo "Run 2"
bash src/stream/run_stream.sh src/stream/stream.10M
echo "Run 3"
bash src/stream/run_stream.sh src/stream/stream.10M
echo "Run 4"
bash src/stream/run_stream.sh src/stream/stream.10M
echo "Run 5"
bash src/stream/run_stream.sh src/stream/stream.10M
echo "Running STREAM_ARRAY_SIZE = 100M" 
echo "Run 1"
bash src/stream/run_stream.sh src/stream/stream.100M
echo "Run 2"
bash src/stream/run_stream.sh src/stream/stream.100M
echo "Run 3"
bash src/stream/run_stream.sh src/stream/stream.100M
echo "Run 4"
bash src/stream/run_stream.sh src/stream/stream.100M
echo "Run 5"
bash src/stream/run_stream.sh src/stream/stream.100M


echo "WARNING!!!" 
echo "  After inspecting the containts of ${WORKING_PATH}/benchmarking_out, you may "
echo "  consider deleteing it"
