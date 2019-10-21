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
import numpy as np
import random
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
            "USAGE : ./src/numpy_attempt.py \n"
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
    elif(nArg != 1):
        print_help(1)
    ### Variables ###
    matrixSize = 5000     # outer dim of mats to run matrix_multiply on
    N = 50       # shortened num of trials to test,get stdev and mean

    random.seed(42)
    np.random.seed(42)
    Ax = matrixSize
    Ay = 10000
    Bx = 10000
    By = matrixSize
    A=np.random.rand(Ax,Ay)
    B=np.random.rand(Bx,By)
    print("Matrix Size = [{} {}] x [{} {}] = multiplied {} times ".format(Ax,Ay,Bx,By,N))
    npStartTime = time.time()
    for i in range(N):
        AB = np.dot(A,B)
    print("Run time : {:.4f} s".format((time.time() - npStartTime)))


if __name__ == "__main__":
    main()
