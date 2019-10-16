# Author : Ali Snedden
# Date   : 5/15/19
# License: MIT
# Purpose: 
#
# Notes :
#
# Questions:
#
# References :  
#
import time
import sys
import os
import glob
import numpy as np
import subprocess
import shutil
import bisect
from error import exit_with_error
from functions import parse_run_time
from functions import get_core_ids

def print_help(Arg):
    """
    ARGS:
        arg     : exit value
    RETURN:
        N/A
    DESCRIPTION:
        Print Help. Exit with value arg
    DEBUG:
        1. Tested, it worked
    FUTURE:
    """
    sys.stdout.write(
            "USAGE : ./src/driver.py [options] /abs/path/to/workspace  /abs/path/to/ref\n"
            "   [options] = 'all' : Builds data, runs all tests\n"
            "             = 'build_mat_mult_data' : Builds matrix_multiply data \n"
            "   NOTES : \n"
            )
    sys.exit(Arg)


def main():
    """
    ARGS:

    RETURN:
    DESCRIPTION:
    NOTES:
    DEBUG:
    FUTURE:
    """
    ### Check Python version and CL args ###
    if(sys.version_info[0] != 3):
        exit_with_error("ERROR!!! Runs with python3, NOT python-{}\n\n".format(
                        sys.version_info[0]))
    nArg = len(sys.argv)
    if(nArg == 2 and (sys.argv[1][0:3] == "--h" or sys.argv[1][0:2] == "-h")):
        print_help(0)
    elif(nArg != 3):
        print_help(1)
    startTime = time.time()
    print("Start Time : {}".format(time.strftime("%a, %d %b %Y %H:%M:%S ",
                                   time.localtime())))
    print("Logging run output to driver.log\n\n")
    ### Variables ###
    options    = sys.argv[1]
    workPath   = sys.argv[2]           # Path where all the output/work will be saved.
    ompNumThreadsL = [1,2,5,20]        # Cores used in OMP tasks
    matrixSizeL = [10000]     # outer dim of mats to run matrix_multiply on
    #matrixSizeL = [2000,3000,5000]     # outer dim of mats to run matrix_multiply on
    nTrials     = 3                     # number of trials to test,get stdev and mean
    shortNTrials= 1                # shortened num of trials to test,get stdev and mean
    # Create work path dir if doesn't exist
    if(not os.path.isdir(workPath)):
        os.mkdir(workPath)

    if(options != 'all' and options != 'build_mat_mult_data'):
        exit_with_error("ERROR!!! {} is invalid option\n".format(options))


    ######## Run Tests ########
    if(options == 'all' or options == 'build_mat_mult_data'):
        nThread = 1
        #print("Building data for matrix_multiply (time to run is for numpy's matrix mult.: ")
        #print("--------------------------------------------------------")
        #print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
        #      "stdev"))
        #print("--------------------------------------------------------")
        ### Create directory structure in data
        outDirPrefix = "{}/data/matrix".format(workPath)
        if(not os.path.isdir(outDirPrefix)):
            os.mkdir(outDirPrefix)

        for size in matrixSizeL:
            outDir = "{}/{}".format(outDirPrefix,size)
            if(not os.path.isdir(outDir)):
                os.mkdir(outDir)
            runTimeV = np.zeros([shortNTrials])
            for tIdx in range(shortNTrials):   ### change to shortNTrials
                cmd =  ("python3 src/matrix/matrix_generator.py {} 10000 "
                             "10000 {} {}".format(size, size, outDir))
                output = "{}\n".format(cmd)
                output = output + subprocess.getoutput(cmd)
                runTime = parse_run_time(output,workPath) # Run time
                runTimeV[tIdx]= runTime
            print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
            print("--------------------------------------------------------")




if __name__ == "__main__":
    main()
