# Author : Ali Snedden
# Date   : 5/15/19
# Purpose:
#
# Notes :
#
# Questions:
#
import numpy as np
import sys
import subprocess
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



def get_core_ids(NumThreads=None):
    """
    ARGS:
        NumThreads : number of threads
    RETURN:
        List of core IDs
    DESCRIPTION:
        Given a number of threads, return core ids that minimize the number of NUMA
        nodes used and minimize the number of sockets used.
    NOTE: 
        1. I'm making the assumption that the cores are listed sequentially 
           are in the same socket / numa node. This all breaks down if 
           something crazy like this occurs:
        
            coreId | NUMA Node | socket
            ---------------------------
            0      | 0         | 0
            1      | 1         | 0
            2      | 0         | 0
            3      | 1         | 0
            4      | 1         | 0
            5      | 0         | 0
            6      | 0         | 0

          Instead of 

            coreId | NUMA Node | socket
            ---------------------------
            0      | 0         | 0
            1      | 0         | 0
            2      | 0         | 0
            3      | 0         | 0
            4      | 1         | 0
            5      | 1         | 0
            6      | 1         | 0
    DEBUG:
        1. Added sanity checks at the end that seem to work
        2. Tested [1, 2, 5, 6, 7, 10, 12, 15, 20] and returns 
           correct cores as expected
    FUTURE:
    """
    cmd="grep -P 'physical id' /proc/cpuinfo | cut -d: -f2 | tr -d ' '"
    socketL = subprocess.getoutput(cmd).split()
    socketL = [int(socket) for socket in socketL]
    cmd="grep -P 'processor' /proc/cpuinfo | cut -d: -f2 | tr -d ' '"
    coreL = subprocess.getoutput(cmd).split()
    coreL = [int(core) for core in coreL]
    cmd="lscpu | grep 'NUMA node[^(]' | awk '{print $2 \" \" $4}'"
    numaL = subprocess.getoutput(cmd).split('\n')
    cmd="lscpu | grep 'Core(s) per socket' | awk '{print $4}'"
    nCorePerSocket = int(subprocess.getoutput(cmd))

    coreNumaD=dict()
    coreSockD=dict()
    for node in numaL:
        [node,coreRange] = node.split()
        coreRange=coreRange.split('-')
        nCorePerNuma = int(coreRange[1]) -  int(coreRange[0]) + 1 # Assuming same cores

        for coreId in range(int(coreRange[0]), int(coreRange[1]) + 1):
            # Get core's NUMA node
            coreNumaD[coreId] = node

            # Get core's Socket
            for idx in range(len(coreL)):
                if(coreId == coreL[idx]):
                    coreSockD[coreId] = socketL[idx]

    ## If asking for more threads then cores, just max to number of cores
    if(len(coreL) < NumThreads):
        coreId2returnL = coreL
        return(coreId2returnL)
    else:
        coreId2returnL = np.arange(min(coreL), NumThreads)

    # Error Check
    if(len(coreId2returnL) != NumThreads):
        exit_with_error("ERROR!!! coreId2returnL {} != NumThreads "
                        "{}".format(coreId2returnL,NumThreads))
    # Check that 
    retSocketL = []
    retNumaL   = []
    for core in coreId2returnL:
        # Numa
        if(coreNumaD[core] not in retNumaL): 
            retNumaL.append(coreNumaD[core])
        # Socket
        if(coreSockD[core] not in retSocketL): 
            retSocketL.append(coreSockD[core])
    # Check for minimum number of NUMA nodes
    if(int(len(coreId2returnL) / (nCorePerNuma*1.001) + 1) != len(retNumaL)):
        exit_with_error("ERROR!!! len(coreId2returnL) {} / nCorePerNuma {}"
                        " + 1 != len(retNumL) {}\n".format(len(coreId2returnL),
                        nCorePerNuma, len(retNumaL)))
    # Check for minimum number of sockets
    if(int(len(coreId2returnL) / (nCorePerSocket*1.001) + 1) != len(retSocketL)):
        exit_with_error("ERROR!!! len(coreId2returnL) {} / nCorePerSocket {}"
                        " + 1 != len(retSocketL) {}\n".format(len(coreId2returnL),
                        nCorePerSocket, len(retSocketL)))
    
    return(coreId2returnL)
