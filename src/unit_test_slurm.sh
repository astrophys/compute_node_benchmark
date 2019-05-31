#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --output=%x_%j.out
# Author : Ali Snedden
# Date : 5/30/19
# License: MIT
# Purpose : 
#   This batch script helps you unit test the singularity image using
#   the SLURM scheduler. This makes it possible to actually test the 
#   components of the container in a fashion that takes much less than 24h
#
# Run : 
#       sbatch submit.slrm [unit_test_option]
#
#
# Comprehensibly run :
#       for unit in `echo build_mat_mult_data mat_mult_cache_opt mat_mult_non_cache_opt build_rnaseq_data align_rnaseq_tophat align_rnaseq_hisat cufflinks_assemble cuffmerge cuffcompare cuffquant cuffnorm cuffdiff kelvin stream`; do sbatch --mail-user=your@email.com --mail-type=FAIL  --job-name=${unit} unit_test_slurm.sh $unit; done
UNIT_TEST=$1


## Test all
module load singularity/2.5.1
if [[ ! -z "${UNIT_TEST}" ]]; then
    if [ -d "/tmp/nodeTest" ]; then
        rm -r /tmp/nodeTest
    fi
    cp -rp /gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking/nodeTest /tmp
    PREV_OUTPUT=/tmp/nodeTest


    ## For most of these unit tests, previous output must have been generated
    ## so we'll dump the current output to the same dir as PREV_OUTPUT
    if [[ -z "${PREV_OUTPUT}" ]]; then
        echo "ERROR!!! UNIT_TEST set but PREV_OUTPUT is unset" >&2
        exit 1
    fi

    ### Test if ok
    if [ ${UNIT_TEST} != "build_mat_mult_data"] && \
       [ ${UNIT_TEST} != "mat_mult_cache_opt" ] && \
       [ ${UNIT_TEST} != "mat_mult_non_cache_opt" ] && \
       [ ${UNIT_TEST} != "build_rnaseq_data" ] && \
       [ ${UNIT_TEST} != "align_rnaseq_tophat" ] && \
       [ ${UNIT_TEST} != "align_rnaseq_hisat" ] && \
       [ ${UNIT_TEST} != "cufflinks_assemble" ] && \
       [ ${UNIT_TEST} != "cuffmerge" ] && \
       [ ${UNIT_TEST} != "cuffcompare" ] && \
       [ ${UNIT_TEST} != "cuffquant" ] && \
       [ ${UNIT_TEST} != "cuffnorm" ] && \
       [ ${UNIT_TEST} != "cuffdiff" ] && \
       [ ${UNIT_TEST} != "kelvin" ] && \
       [ ${UNIT_TEST} != "stream" ]; then
        echo "ERROR!!! UNIT_TEST=${UNIT_TEST} is an invalid option"
        exit 1
    else
        if [ ${UNIT_TEST} == "build_mat_mult_data" ]; then
            rm -r ${PREV_OUTPUT}/data/matrix
        elif [ ${UNIT_TEST} == "mat_mult_cache_opt" ]; then
            rm -r ${PREV_OUTPUT}/output/matrix_cache_opt
        elif [ ${UNIT_TEST} == "mat_mult_non_cache_opt" ]; then
            rm -r ${PREV_OUTPUT}/output/matrix_non_cache_opt
        elif [ ${UNIT_TEST} == "build_rnaseq_data" ]; then
            rm -r ${PREV_OUTPUT}/data/rnaseq
        elif [ ${UNIT_TEST} == "align_rnaseq_tophat" ]; then
            rm -r ${PREV_OUTPUT}/output/rnaseq/tophat
        elif [ ${UNIT_TEST} == "align_rnaseq_hisat" ]; then
            rm -r ${PREV_OUTPUT}/output/rnaseq/hisat
        elif [ ${UNIT_TEST} == "cufflinks_assemble" ]; then
            rm -r ${PREV_OUTPUT}/output/rnaseq/cufflinks
        elif [ ${UNIT_TEST} == "cuffmerge" ]; then
            rm -r ${PREV_OUTPUT}/output/rnaseq/cuffmerge
        elif [ ${UNIT_TEST} == "cuffcompare" ]; then
            rm -r ${PREV_OUTPUT}/output/rnaseq/cuffcompare
        elif [ ${UNIT_TEST} == "cuffquant" ]; then
            rm -r ${PREV_OUTPUT}/output/rnaseq/cuffquant
        elif [ ${UNIT_TEST} == "cuffnorm" ]; then
            rm -r ${PREV_OUTPUT}/output/rnaseq/cuffnorm
        elif [ ${UNIT_TEST} == "cuffdiff" ]; then
            rm -r ${PREV_OUTPUT}/output/rnaseq/cuffdiff
        elif [ ${UNIT_TEST} == "kelvin" ]; then
            rm -r ${PREV_OUTPUT}/output/kelvin
        elif [ ${UNIT_TEST} == "stream" ]; then
            echo "Nothing done"
        else
            echo "ERROR!!! UNIT_TEST = ${UNIT_TEST} is an invalid option"
            exit 1
        fi
    fi
fi

## Unit Test
export SINGULARITYENV_UNIT_TEST=${UNIT_TEST}
export SINGULARITYENV_PREV_OUTPUT=${PREV_OUTPUT}
singularity run -H /home/gdhpcgroup/aps003 test.simg   ## Runs prebuilt tests


