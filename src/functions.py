# Author : Ali Snedden
# Date   : 5/15/19
# Purpose:
#
# Notes :
#
# Questions:
#
import sys
from error import exit_with_error


def parse_run_time(String, OutDir):
    """
    ARGS:
        String : A string containing the output of run
        OutDir : This matters b/c we need to point to a place 
                 where our Singularity container can actually write
    RETURN:
        Run Time from the output in String
    DESCRIPTION:
        Parses String for "Run time : " . Captures time and returns
        it as a float. Tries to catch errors. 

        WARNING : This function is delicately sensitive to the exact
                  format of the "Run time : " string. It is not robust
                  from this aspect. It also excepts the output from Linux time
                  "real    "
                  
    DEBUG:
        1. Spot checked that it can parse the run time
        2. Checked that both error check conditionals work
    FUTURE:
    """
    stringL = String.split("\n")
    time = None
    # Report log of run
    driverLogF = open("{}/driver.log".format(OutDir), "a")
    driverLogF.write("\n\n\n\n##################################################\n")
    driverLogF.write(String)
    driverLogF.close()

    # Parse
    for string in stringL:
        if("ERROR" in string or "command not found" in string or 
           "Errno" in string
        ):
            exit_with_error("{}\n".format(string))
        elif("Run time :" in string):
            time = float(string.split(" ")[3])
            break
        elif("real\t" in string):
            string = string.split()[1]
            string = string.split('m')
            minutes= float(string[0])
            seconds= float(string[1].split('s')[0])
            time = minutes*60 + seconds
            break
    if(time is None):
        exit_with_error("ERROR!!! Incorrectly formatted 'Run time : '\n".format(string))

    return(time)



