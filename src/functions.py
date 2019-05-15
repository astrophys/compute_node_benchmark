# Author : Ali Snedden
# Date   : 5/15/19
# Purpose:
#
# Notes :
#
# Questions:
#
import time
import sys
from error import exit_with_error


def parse_run_time(String):
    """
    ARGS:
        String : A string containing the output of run
    RETURN:
        Run Time from the output in String
    DESCRIPTION:
        Parses String for "Run time : " . Captures time and returns
        it as a float. Tries to catch errors. 

        WARNING : This function is delicately sensitive to the exact
                  format of the "Run time : " string. It is not robust
                  from this aspect.
    DEBUG:
        1. Spot checked that it can parse the run time
        2. Checked that both error check conditionals work
    FUTURE:
    """
    stringL = String.split("\n")
    time = None
    for string in stringL:
        if("ERROR" in string):
            exit_with_error("{}\n".format(string))
        elif("Run time :" in string):
            time = float(string.split(" ")[3])
            break
    if(time is None):
        exit_with_error("ERROR!!! Incorrectly formatted 'Run time : '\n".format(string))

    return(time)



