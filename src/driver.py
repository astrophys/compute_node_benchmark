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
import numpy as np
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
            "USAGE : ./src/driver.py [options]\n"
            "   [options] = 'all' : Builds data, runs all tests\n"
            "             = 'build_mat_mult_data'    : Builds matrix_multiply data \n"
            "             = 'mat_mult_cache_opt'     : Run matrix_multiply_cache_opt tests\n"
            "             = 'mat_mult_non_cache_opt' : Run matrix_multiply_cache_opt tests\n"
            "             = 'build_rnaseq_data'      : Creates single end RNA-Seq data\n"
            "             = 'local_memory_access'    : grep's a large file in temp\n"
            "   NOTE : Only one option can be passed at a time.\n"
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
    elif(nArg != 2):
        print_help(1)
    startTime = time.time()
    print("Start Time : {}".format(time.strftime("%a, %d %b %Y %H:%M:%S ",
                                   time.localtime())))
    print("Logging run output to driver.log\n\n")
    ### Variables ###
    options     = sys.argv[1]
    ompNumThreadsL = [1,2,5,7,10,15,20]  ## Cores used in OMP tasks
    matrixSizeL = [2000,3000,5000]       ## outer dim of mats to run matrix_multiply on
    #rnaSeqSizeL = [10**4,10**5,10**6]   
    rnaSeqSizeL = [2*10**4,10**5]   
    nTrials     = 3                      ## number of trials to test,get stdev and mean

    if(options != 'all' and options != 'mat_mult_cache_opt' and 
       options != 'mat_mult_non_cache_opt' and options != 'local_memory_access' and
       options != 'build_mat_mult_data' and options != 'build_rnaseq_data' 
    ):
        exit_with_error("ERROR!!! {} is invalid option\n")


    ######## Run Tests ########
    if(options == 'all' or options == 'build_mat_mult_data'):
        nThread = 1
        print("building data for matrix_multiply : ")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        for size in matrixSizeL:
            outDir = "data/{}".format(size)
            runTimeV = np.zeros([nTrials])
            for tIdx in range(nTrials):
                cmd =  "python3 src/matrix/matrix_generator.py {} 10000 10000 {} {}".format(
                       size, size, outDir)
                output = subprocess.getoutput(cmd)
                runTime = parse_run_time(output) # Run time
                runTimeV[tIdx]= runTime
            print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
            print("--------------------------------------------------------")


    if(options == 'all' or options == 'mat_mult_cache_opt'):
        print("matrix_multiply (cache optimized using OpenMP) : ")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        for size in matrixSizeL:
            for nThread in ompNumThreadsL:
                runTimeV = np.zeros([nTrials])
                #nThread = 10
                #size=2000
                for tIdx in range(nTrials):
                    cmd =  ("export OMP_NUM_THREADS={}; ./src/matrix/matrix_multiply_cache_opt "
                            "data/{}/A.txt data/{}/B.txt  "
                             "data/{}/output".format(nThread,size,size,size))
                    output = subprocess.getoutput(cmd)
                    runTime = parse_run_time(output) # Run time
                    runTimeV[tIdx]= runTime
                print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
            print("--------------------------------------------------------")


    if(options == 'all' or options == 'mat_mult_non_cache_opt'):
        print("matrix_multiply (non-cache optimized using OpenMP) : ")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        for size in matrixSizeL:
            for nThread in ompNumThreadsL:
                runTimeV = np.zeros([nTrials])
                #nThread = 10
                #size=2000
                for tIdx in range(nTrials):
                    cmd =  ("export OMP_NUM_THREADS={}; "
                            "./src/matrix/matrix_multiply_non_cache_opt "
                            "data/{}/A.txt data/{}/B.txt  "
                             "data/{}/output".format(nThread,size,size,size))
                    output = subprocess.getoutput(cmd)
                    runTime = parse_run_time(output) # Run time
                    runTimeV[tIdx]= runTime
                print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
            print("--------------------------------------------------------")


    if(options == 'all' or options == 'build_rnaseq_data'):
        print("Building RNA-Seq Data sets : ")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        nThread = 1
        nSamp   = 3
        outDir = "~/Code/Singularity/benchmarking/data/rnaseq/fastq/"
        gtf="~/Code/Singularity/benchmarking/src/simulate_fastq_data/data/chr1_short.gtf"
        genome ="~/Code/Singularity/benchmarking/src/simulate_fastq_data/data/chr1_short.fa" 
        configL=["config/config_wt_chr1.txt", "config/config_treat_chr1.txt"]
        #gtf="/reference/homo_sapiens/GRCh38/ensembl/release-83/Annotation/Genes/gtf/Homo_sapiens.GRCh38.83.gtf"
        #genome ="/reference/homo_sapiens/GRCh38/ensembl/release-83/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        #configL=["config/config_wt.txt", "config/config_treat.txt"]
        ## Loop
        for size in rnaSeqSizeL:
            runTimeV = np.zeros([nSamp*len(configL)])
            tIdx = 0
            for config in configL:
                for samp in range(nSamp):
                    ## Set output files
                    if("treat" in config):
                        outFile = "{}/treat_{}_{}".format(outDir,size,samp)
                    elif("wt" in config):
                        outFile = "{}/wt_{}_{}".format(outDir,size,samp)
                    else:
                        exit_with_error("ERROR!!! No correct configuration file found!\n")
                    cmd =  ("export OMP_NUM_THREADS={}; "
                       "python3 src/simulate_fastq_data/simulate_fastq.py "
                       "{} {} {} {} {} single"
                       "".format(nThread, gtf, genome, config, size, outFile))
                    output = subprocess.getoutput(cmd)
                    runTime = parse_run_time(output) # Run time
                    runTimeV[tIdx]= runTime
                    tIdx = tIdx + 1
            print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                  np.mean(runTimeV), np.std(runTimeV)))
            print("--------------------------------------------------------")





    print("Run Time : {:.4f} h".format((time.time() - startTime)/3600.0))


if __name__ == "__main__":
    main()
