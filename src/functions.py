# Author : Ali Snedden
# Date   : 5/15/19
# License: MIT
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
    driverLogF.flush()
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
        1. Previously, I made the assumption that the cores are listed sequentially 
           are in the same socket / numa node. This all breaks down if 
           something crazy happens (like cores on alternating NUMA nodes).  I believe
           that I have fixed for this.
        2. It is ok that I don't sort coreId2returnL b/c taskset doesn't care
           about the order.
        
    DEBUG:
        1. Added sanity checks at the end that seem to work
        2. Tested [1, 2, 5, 6, 7, 10, 12, 15, 20] and returns 
           correct cores as expected
        3. Added output from lscpu that DELL provided, changed cmd appropriately
           to emulate how lscpu returns values.

           Added loop to handle cores being on alternating socket/numa
           --> Tested interleave_off.txt, [1,2,10,20,25,40,50], returns as expected
           --> Tested interleave_on.txt,  [1,2,10,20,25,40,50], returns as expected
           --> Retested on Baker, but himem and gp nodes: [1,3,5,6,12,24,48,100],
               returns as expected.
    FUTURE:
    """
    # Create several lists where the elements of the lists correspond to the individual
    # cores and their properties
    #       E.g. socketL[1] corresponds the the socket core (with coreID = coreL[1])
    #cmd="grep -P 'physical id' /proc/cpuinfo | cut -d: -f2 | tr -d ' '"
    cmd="lscpu -e | grep -v CORE | awk '{print $3}'"
    socketL = subprocess.getoutput(cmd).split()
    socketL = [int(socket) for socket in socketL]   # socket per list
    #cmd="grep -P 'processor' /proc/cpuinfo | cut -d: -f2 | tr -d ' '"
    cmd="lscpu -e | grep -v CORE | awk '{print $4}'"
    coreL = subprocess.getoutput(cmd).split()
    coreL = [int(core) for core in coreL]
    #cmd="lscpu -e | grep 'NUMA node[^(]' | awk '{print $2 \" \" $4}'"
    #numaL = subprocess.getoutput(cmd).split('\n')
    cmd="lscpu -e | grep -v CORE | awk '{print $2}'"
    numaL = subprocess.getoutput(cmd).split()
    numaL = [int(numa) for numa in numaL]       

    cmd="lscpu | grep 'Core(s) per socket' | awk '{{print $4}}'"
    nCorePerSocket = int(subprocess.getoutput(cmd))

    cmd="lscpu -e | grep -v CORE | awk '{{print $2}}' | sort -n  | uniq -c | awk '{{print $1}}' | head -n 1"
    nCorePerNuma = int(subprocess.getoutput(cmd).split()[0])

    ####### BEGIN debugging DELL's lscpu ########
    ## Create several lists where the elements of the lists correspond to the individual
    ## cores and their properties
    ##       E.g. socketL[1] corresponds the the socket core (with coreID = coreL[1])
    ## Meant to be called when debugging from within src/
    ##
    ##
    ##fin="interleave_off.txt"
    #fin="interleave_on.txt"
    #cmd="grep -v CORE ../config/dell/{} | awk '{{print $3}}'".format(fin)
    #socketL = subprocess.getoutput(cmd).split()
    #socketL = [int(socket) for socket in socketL]   # socket per list
    ##cmd="grep -P 'processor' /proc/cpuinfo | cut -d: -f2 | tr -d ' '"
    #cmd="grep -v CORE ../config/dell/{} | awk '{{print $4}}'".format(fin)
    #coreL = subprocess.getoutput(cmd).split()
    #coreL = [int(core) for core in coreL]
    ##cmd="lscpu -e | grep 'NUMA node[^(]' | awk '{print $2 \" \" $4}'"
    ##numaL = subprocess.getoutput(cmd).split('\n')
    #cmd="grep -v CORE ../config/dell/{} | awk '{{print $2}}'".format(fin)
    #numaL = subprocess.getoutput(cmd).split()
    #numaL = [int(numa) for numa in numaL]       ### This changed!
    #cmd="grep 'Core(s) per socket' ../config/dell/lscpu_{} | awk '{{print $4}}'".format(fin)
    #nCorePerSocket = int(subprocess.getoutput(cmd))
    #cmd="grep -v CORE ../config/dell/{} | awk '{{print $2}}' | sort -n  | uniq -c | awk '{{print $1}}' | head -n 1".format(fin)
    #nCorePerNuma = int(subprocess.getoutput(cmd).split()[0])
    ####### END debugging DELL's lscpu ########

    coreNumaD=dict()
    coreSockD=dict()
    ## Create dicts 
    for idx in range(len(coreL)):
        # Get core's NUMA node
        coreId = coreL[idx]
        coreNumaD[coreId] = numaL[idx]
        # Get core's Socket
        coreSockD[coreId] = socketL[idx]

    #import pdb; pdb.set_trace()
    ## If asking for more threads then cores, just max to number of cores
    if(len(coreL) < NumThreads):
        coreId2returnL = coreL
        string=",".join(str(i) for i in coreId2returnL)
        return(string)
    else:
        n = 0              ## Number of cores found
        numaNode = 0       ## Current NUMA node working on
        socket   = 0
        coreId2returnL =[] ## Cores found
        while(n<NumThreads):
            for idx in range(len(coreL)):
                if(numaL[idx] == numaNode and socketL[idx] == socket):
                    coreId2returnL.append(coreL[idx])
                    n = n + 1
                if(n == NumThreads):
                    break

            if(n >= nCorePerNuma*(numaNode + 1)):
                numaNode = numaNode + 1
            # Iterate socket as run out of space
            if(n >= nCorePerSocket*(socket + 1)):
                socket = socket + 1
            
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
        exit_with_error("ERROR!!! len(coreId2returnL) ({}) / nCorePerNuma {}"
                        " + 1 != len(retNumL) ({})\n".format(len(coreId2returnL),
                        nCorePerNuma, len(retNumaL)))
    # Check for minimum number of sockets
    if(int(len(coreId2returnL) / (nCorePerSocket*1.001) + 1) != len(retSocketL)):
        exit_with_error("ERROR!!! len(coreId2returnL) {} / nCorePerSocket {}"
                        " + 1 != len(retSocketL) {}\n".format(len(coreId2returnL),
                        nCorePerSocket, len(retSocketL)))
    string=",".join(str(i) for i in coreId2returnL)
    return(string)
