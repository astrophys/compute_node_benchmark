# Author : Ali Snedden
# Date   : 5/15/19
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
import subprocess
from error import exit_with_error
from functions import parse_run_time

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
            "USAGE : ./src/driver.py "
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
    startTime = time.time()
    print("Start Time : {}".format(time.strftime("%a, %d %b %Y %H:%M:%S ",
                                   time.localtime())))

    nCoresL = [1,2,5,7,10,15,20] ## Cores used in OMP tasks
    matrixSizeL = [2000,3000,5000,10000]
    ## Run Tests ##
    #for nCores in omp
    cmd =  "export OMP_NUM_THREADS=10; ./src/matrix/matrix_multiply data/2000/A.txt data/2000/B.txt  data/2000/output" 
    output = subprocess.getoutput(cmd)
    runTime = parse_run_time(output) # Run time
    print(runTime)

    print("Run Time : {:.4f} h".format((time.time() - startTime)/3600.0))


if __name__ == "__main__":
    main()
