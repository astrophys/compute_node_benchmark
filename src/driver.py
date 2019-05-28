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
import os
import glob
import numpy as np
import subprocess
import shutil
import bisect
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
            "USAGE : ./src/driver.py [options] /abs/path/to/workspace  /abs/path/to/ref\n"
            "   [options] = 'all' : Builds data, runs all tests\n"
            "             = 'build_mat_mult_data' : Builds matrix_multiply data \n"
            "             = 'mat_mult_cache_opt'  : Run matrix_multiply_cache_opt tests\n"
            "             = 'mat_mult_non_cache_opt' : Run matrix_multiply_cache_opt tests\n"
            "             = 'build_rnaseq_data'   : Creates single end RNA-Seq data\n"
            "             = 'align_rnaseq_tophat' : Align samples in data/rnaseq w/ tophat\n"
            "             = 'align_rnaseq_hisat'  : Align samples in data/rnaseq w/ hisat\n"
            "             = 'cufflinks_assemble'  : Must have run tophat. Assembles transcriptome\n"
            "             = 'cuffmerge'           : Must have run tophat and cufflinks\n"
            "             = 'cuffcompare'         : Must have run tophat,cufflinks\n"
            "             = 'cuffquant'           : Must have run tophat,cufflinks,cuffmerge\n"
            "             = 'cuffnorm'            : Must have run tophat,cufflinks,"
            "             = 'cuffdiff'            : Must have run tophat,cufflinks,"
                                                    "cuffmerge and cuffquant\n"
            "             = 'kelvin'              : Runs kelvin (a statistical genetics software) \n"
            "             = 'local_memory_access' : grep's a large file in temp\n"
            "   /abs/path/to/workspace    : Path where all the output/data gets saved.\n"
            "   /abs/path/to/ref          : Path where bowtie/hisat indices and ref fasta/gtf"
            "                               are stored\n\n"
            "   NOTES : \n"
            "       1. Only one option can be passed at a time.\n"
            "       2. It is assumed that all exectuables (i.e. tophat2, bowtie2, etc.)\n"
            "          are located in you shell PATH\n"
            "       3. Location / Names of references _must_ be :\n"
            "          a) Bowtie2 indices: /refPath/Bowtie2Index/Homo_sapiens.GRC38 \n"
            "          b) Hisat2 indices : /refPath/HisatIndex/genome \n"
            "          c) Genome Fasta   : /refPath/Homo_sapiens.GRCh38.dna.primary_assembly.fa:\n"
            "          d) GTF file       : /refPath/Homo_sapiens.GRCh38.83.gtf\n"
            "          e) Short Chr1 Gtf : /refPath/chr1_short.gtf\n"
            "          f) Short Chr1 fasta : /refPath/chr1_short.fa\n"
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
    elif(nArg != 4):
        print_help(1)
    startTime = time.time()
    print("Start Time : {}".format(time.strftime("%a, %d %b %Y %H:%M:%S ",
                                   time.localtime())))
    print("Logging run output to driver.log\n\n")
    ### Variables ###
    options    = sys.argv[1]
    workPath   = sys.argv[2]           # Path where all the output/work will be saved.
    refPath    = sys.argv[3]           # Path where all the ref data and indices are located
    ompNumThreadsL = [1,2,5,7,10,15,20] # Cores used in OMP tasks
    matrixSizeL = [2000,3000,5000]      # outer dim of mats to run matrix_multiply on
    #matrixSizeL = [2000,3000]      # outer dim of mats to run matrix_multiply on
    #rnaSeqSizeL = [10**4,10**5,10**6]   
    rnaSeqSizeL = [2*10**4,10**5]
    nTrials     = 3                     # number of trials to test,get stdev and mean
    # Create work path dir if doesn't exist
    if(not os.path.isdir(workPath)):
        os.mkdir(workPath)
    # Get operating system and list of cores (linux only)
    curOS = sys.platform
    if(curOS == 'darwin'):
        curOS = 'osx'                   # Rename for my own selfish readability
    elif(curOS == 'linux'):
        cmd = "grep -P 'processor[\t ]' /proc/cpuinfo | cut -d: -f2 | tr -d ' '"
        coreIDL = subprocess.getoutput(cmd)
        coreIDL = [int(idx) for idx in coreIDL.split()]
    else:
        exit_with_error("ERROR!! {} is an unsupported operating system".format(curOS))


    ## In Linux singularity container add cores per socket and total cores to ompNumThreadsL
    if(shutil.which('lscpu') != None):
        cmd="lscpu | grep 'Core(s) per socket:' | awk '{print $4}'"
        coresPerSocket = int(subprocess.getoutput(cmd))
        cmd="lscpu  | grep '^CPU(s):' | awk '{print $2}'"
        totalCores = int(subprocess.getoutput(cmd))
        bisect.insort_left(ompNumThreadsL, coresPerSocket)
        bisect.insort_left(ompNumThreadsL, totalCores)
        ompNumThreadsL=list(sorted(set(ompNumThreadsL)))
        print("Cores per socket : {}".format(coresPerSocket))
        print("Total Cores : {}".format(totalCores))
        print("Cores tested : {}".format(ompNumThreadsL))

    if(options != 'all' and options != 'build_mat_mult_data' and 
       options != 'mat_mult_non_cache_opt' and options != 'local_memory_access' and
       options != 'mat_mult_cache_opt' and options != 'build_rnaseq_data' and
       options != 'align_rnaseq_tophat' and options != 'align_rnaseq_hisat' and
       options != 'cufflinks_assemble'  and options != 'cuffmerge' and
       options != 'cuffcompare' and options != 'cuffquant' and
       options != 'cuffnorm' and options != 'cuffdiff' and options != 'kelvin'
    ):
        exit_with_error("ERROR!!! {} is invalid option\n".format(options))


    ######## Run Tests ########
    if(options == 'all' or options == 'build_mat_mult_data'):
        nThread = 1
        print("Building data for matrix_multiply (time to run is for numpy's matrix mult.: ")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        ### Create directory structure in data
        outDirPrefix = "{}/data/matrix".format(workPath)
        if(not os.path.isdir(outDirPrefix)):
            os.mkdir(outDirPrefix)

        for size in matrixSizeL:
            outDir = "{}/{}".format(outDirPrefix,size)
            if(not os.path.isdir(outDir)):
                os.mkdir(outDir)
            runTimeV = np.zeros([nTrials])
            for tIdx in range(nTrials):   ### change to nTrials
                cmd =  "python3 src/matrix/matrix_generator.py {} 10000 10000 {} {}".format(
                       size, size, outDir)
                output = "{}\n".format(cmd)
                output = output + subprocess.getoutput(cmd)
                runTime = parse_run_time(output,workPath) # Run time
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

        ### Create directory structure in output
        outDirPrefix = "{}/output/matrix_cache_opt".format(workPath)
        if(not os.path.isdir(outDirPrefix)):
            os.mkdir(outDirPrefix)

        for size in matrixSizeL:
            outDir = "{}/{}".format(outDirPrefix,size)
            if(not os.path.isdir(outDir)):
                os.mkdir(outDir)

            for nThread in ompNumThreadsL:
                runTimeV = np.zeros([nTrials])
                #nThread = 10
                #size=2000
                for tIdx in range(nTrials):
                    cmd =  ("export OMP_NUM_THREADS={}; ./src/matrix/matrix_multiply_cache_opt "
                            "{}/data/matrix/{}/A.txt {}/data/matrix/{}/B.txt  "
                             "{}".format(nThread,workPath,size,workPath,size,outDir))
                    output = "{}\n".format(cmd)
                    output = output + subprocess.getoutput(cmd)
                    runTime = parse_run_time(output,workPath) # Run time
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

        ### Create directory structure in output
        outDirPrefix = "{}/output/matrix_non_cache_opt".format(workPath)
        if(not os.path.isdir(outDirPrefix)):
            os.mkdir(outDirPrefix)

        for size in matrixSizeL:
            outDir = "{}/{}".format(outDirPrefix,size)
            if(not os.path.isdir(outDir)):
                os.mkdir(outDir)

            for nThread in ompNumThreadsL:
                runTimeV = np.zeros([nTrials])
                #nThread = 10
                #size=2000
                for tIdx in range(nTrials):
                    cmd =  ("export OMP_NUM_THREADS={}; "
                            "./src/matrix/matrix_multiply_non_cache_opt "
                            "{}/data/matrix/{}/A.txt {}/data/matrix/{}/B.txt  "
                             "{}".format(nThread,workPath,size,workPath,size,outDir))
                    output = "{}\n".format(cmd)
                    output = output + subprocess.getoutput(cmd)
                    runTime = parse_run_time(output,workPath) # Run time
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
        treatSampL = []
        wtSampL = []
        gtf="{}/chr1_short.gtf".format(refPath)
        genome ="{}/chr1_short.fa".format(refPath)
        configL=["config/config_wt_chr1.txt", "config/config_treat_chr1.txt"]

        # Create output directory structure
        outDir  = "{}/data/rnaseq".format(workPath)
        if(not os.path.isdir(outDir)):
            os.mkdir(outDir)
        outDir  = "{}/fastq/".format(outDir)
        if(not os.path.isdir(outDir)):
            os.mkdir(outDir)
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
                        if(not os.path.isdir("{}/{}".format(outDir,size))):
                            os.mkdir("{}/{}".format(outDir,size))
                        outFile = "{}/{}/treat_{}".format(outDir,size,samp)
                        treatSampL.append(outFile)
                    elif("wt" in config):
                        if(not os.path.isdir("{}/{}".format(outDir,size))):
                            os.mkdir("{}/{}".format(outDir,size))
                        outFile = "{}/{}/wt_{}".format(outDir,size,samp)
                        wtSampL.append(outFile)
                    else:
                        exit_with_error("ERROR!!! No correct config file found!\n")
                    cmd =  ("export OMP_NUM_THREADS={}; "
                       "python3 src/simulate_fastq_data/simulate_fastq.py "
                       "{} {} {} {} {} single"
                       "".format(nThread, gtf, genome, config, size, outFile))
                    output = "{}\n".format(cmd)
                    output = output + subprocess.getoutput(cmd)
                    runTime = parse_run_time(output,workPath) # Run time
                    runTimeV[tIdx]= runTime
                    tIdx = tIdx + 1
            print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                  np.mean(runTimeV), np.std(runTimeV)))
            print("--------------------------------------------------------")


    if(options == 'all' or options == 'align_rnaseq_tophat'):
        print("Aligning RNA-Seq Data sets with tophat : ")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        outDirPref = "{}/output/rnaseq".format(workPath)
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        outDirPref = os.path.abspath("{}/output/rnaseq/tophat".format(workPath)) 
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)

        inDirPref  = os.path.abspath("{}/data/rnaseq/fastq".format(workPath)) 
        if(not os.path.isdir(inDirPref)):
            exit_with_error("ERROR!!! fastq data does not exits. Run build_rnaseq_data option")
        bowtieIdxPath = "{}/Bowtie2Index/Homo_sapiens.GRC38".format(refPath)
        ## Loop
        for size in rnaSeqSizeL:
            sampFileL   = glob.glob("{}/{}/*.fq".format(inDirPref,size))
            if(not os.path.isdir("{}/{}".format(outDirPref,size))):
                os.mkdir("{}/{}".format(outDirPref,size))

            for nThread in ompNumThreadsL:
                runTimeV = np.zeros([len(sampFileL)])
                tIdx = 0
                for samp in sampFileL:
                    sampDir = samp.split("/")[-1].split(".")[0]
                    ## Set output directory
                    outDir = "{}/{}/{}".format(outDirPref,size,sampDir)
                    if(os.path.isfile("/Users/asnedden/.local/virtualenvs/python2.7/bin/python")):
                        # My OSX configuration b/c I use virtualenv
                        python2="source ~/.local/virtualenvs/python2.7/bin/activate;"
                        cmd =  (
                            "{}; time tophat2 -p {} -o {} {} {}"
                            "".format(python2,nThread, outDir,
                            bowtieIdxPath, samp))
                    else:
                    #    # On CentOS, default python is 2.6.6
                    #    python2="/usr/bin/python"
                        cmd =  (
                            "time tophat2 -p {} -o {} {} {}"
                            "".format(nThread, outDir,
                            bowtieIdxPath, samp))
                    
                    output = "{}\n".format(cmd)
                    output = output + subprocess.getoutput(cmd)
                    runTime = parse_run_time(output,workPath) # Run time
                    runTimeV[tIdx]= runTime
                    tIdx = tIdx + 1
                print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
                print("--------------------------------------------------------")


    if(options == 'all' or options == 'align_rnaseq_hisat'):
        print("Aligning RNA-Seq Data sets with hisat : ")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        # Get directory structure
        outDirPref = "{}/output/rnaseq".format(workPath)
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        outDirPref = os.path.abspath("{}/output/rnaseq/hisat".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        inDirPref  = os.path.abspath("{}/data/rnaseq/fastq".format(workPath))   ## prefix
        if(not os.path.isdir(inDirPref)):
            exit_with_error("ERROR!!! fastq data does not exits. Run build_rnaseq_data option")
        hisatIdxPath = "{}/HisatIndex/genome".format(refPath)
        ## Loop
        for size in rnaSeqSizeL:
            sampFileL   = glob.glob("{}/{}/*.fq".format(inDirPref,size))
            if(not os.path.isdir("{}/{}".format(outDirPref,size))):
                os.mkdir("{}/{}".format(outDirPref,size))

            for nThread in ompNumThreadsL:
                runTimeV = np.zeros([len(sampFileL)])
                tIdx = 0
                for samp in sampFileL:
                    sampDir = samp.split("/")[-1].split(".")[0]
                    ## Set output directory
                    outDir = "{}/{}/{}".format(outDirPref,size,sampDir)
                    if(not os.path.isdir(outDir)):
                        os.mkdir(outDir)
                    cmd =  (
                        "time hisat2 -p {} --phred33 -x {} -U {} -S {}/output.sam"
                       "".format(nThread, hisatIdxPath, samp, outDir))
                    output = "{}\n".format(cmd)
                    output = output + subprocess.getoutput(cmd)
                    runTime = parse_run_time(output,workPath) # Run time
                    runTimeV[tIdx]= runTime
                    tIdx = tIdx + 1
                print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
                print("--------------------------------------------------------")


    if(options == 'all' or options == 'cufflinks_assemble'):
        print("Assembling transcriptome using cufflinks: ")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        outDirPref = os.path.abspath("{}/output/rnaseq/".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        outDirPref = os.path.abspath("{}/output/rnaseq/cufflinks".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        inDirPref  = os.path.abspath("{}/output/rnaseq/tophat".format(workPath))   ## prefix
        gtf="{}/Homo_sapiens.GRCh38.83.gtf".format(refPath)
        ## Loop
        for size in rnaSeqSizeL:
            sampFileL   = glob.glob("{}/{}/*/accepted_hits.bam".format(inDirPref,size))
            if(not os.path.isdir("{}/{}".format(outDirPref,size))):
                os.mkdir("{}/{}".format(outDirPref,size))

            for nThread in ompNumThreadsL:
                runTimeV = np.zeros([len(sampFileL)])
                tIdx = 0
                for samp in sampFileL:
                    sampDir = samp.split("/")[-2].split(".")[0]
                    ## Set output directory
                    outDir = "{}/{}/{}".format(outDirPref,size,sampDir)
                    if(not os.path.isdir(outDir)):
                        os.mkdir(outDir)
                    cmd =  (
                        "time cufflinks --num-threads {} -g {} --output-dir {} {}"
                       "".format(nThread, gtf, outDir, samp))
                    output = "{}\n".format(cmd)
                    output = output + subprocess.getoutput(cmd)
                    runTime = parse_run_time(output,workPath) # Run time
                    runTimeV[tIdx]= runTime
                    tIdx = tIdx + 1
                print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
                print("--------------------------------------------------------")


    if(options == 'all' or options == 'cuffmerge'):
        print("Merging assembled transcriptomes using cuffmerge")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        outDirPref = os.path.abspath("{}/output/rnaseq/".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        outDirPref = os.path.abspath("{}/output/rnaseq/cuffmerge".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        inDirPref  = os.path.abspath("{}/output/rnaseq/cufflinks".format(workPath))   ## prefix
        gtf="{}/Homo_sapiens.GRCh38.83.gtf".format(refPath)
        genome="{}/Homo_sapiens.GRCh38.dna.primary_assembly.fa".format(refPath)
        ## Loop
        for size in rnaSeqSizeL:
            sampFileL   = glob.glob("{}/{}/*/transcripts.gtf".format(inDirPref,size))
            outDir = "{}/{}".format(outDirPref,size)

            if(not os.path.isdir(outDir)):
                os.mkdir(outDir)
            assemblyPath = "{}/assemblies.txt".format(outDir)
            if(not os.path.isfile(assemblyPath)):
                assemblyFile = open(assemblyPath, "w+")
                for samp in sampFileL:
                    assemblyFile.write("{}\n".format(samp))

            for nThread in ompNumThreadsL:
                ## Consider adding nTrials here.
                runTimeV = np.zeros([1])
                tIdx = 0
                if(os.path.isfile("/Users/asnedden/.local/virtualenvs/python2.7/bin/python")):
                    # My OSX configuration b/c I use virtualenv
                    python2="source ~/.local/virtualenvs/python2.7/bin/activate;"
                    cmd =  (
                        "{};"
                        "time  cuffmerge --num-threads {} -o {} "
                        "--ref-gtf {} --ref-sequence {} {}"
                        "".format(python2,nThread, outDir, gtf, genome,
                        assemblyPath))
                else:
                #    # On CentOS, default python is 2.6.6
                #    python2="/usr/bin/python"
                    cmd =  (
                        "time  cuffmerge --num-threads {} -o {} "
                        "--ref-gtf {} --ref-sequence {} {}"
                        "".format(nThread, outDir, gtf, genome, assemblyPath))
                output = "{}\n".format(cmd)
                output = output + subprocess.getoutput(cmd)
                runTime = parse_run_time(output,workPath) # Run time
                runTimeV[tIdx]= runTime
                tIdx = tIdx + 1
                print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
                print("--------------------------------------------------------")


    if(options == 'all' or options == 'cuffcompare'):
        print("Comparing cufflinks gtf using cuffcompare")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        # Check and make directory structure
        outDirPref = os.path.abspath("{}/output/rnaseq/".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            exit_with_error("ERROR!!! Expecting {}/output/rnaseq. Must have run tophat "
                            "and cufflinks prior\n".format(workPath))
        outDirPref = os.path.abspath("{}/output/rnaseq/cuffcompare".format(workPath)) 
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        inDirPref  = os.path.abspath("{}/output/rnaseq/cufflinks".format(workPath))   ## prefix
        gtf="{}/Homo_sapiens.GRCh38.83.gtf".format(refPath)
        genome="{}/Homo_sapiens.GRCh38.dna.primary_assembly.fa".format(refPath)
        nThread = 1
        ## Loop
        for size in rnaSeqSizeL:
            sampFileL   = glob.glob("{}/{}/*/transcripts.gtf".format(inDirPref,size))
            outPref = "{}/{}".format(outDirPref,size)

            ## Consider adding nTrials here.
            runTimeV = np.zeros([1])
            tIdx = 0
            cmd =  (
                    "time cuffcompare -o {} -r {} -R -C -V {}"
                    "".format(outPref, gtf, " ".join(sampFileL)))
            output = "{}\n".format(cmd)
            output = output + subprocess.getoutput(cmd)
            runTime = parse_run_time(output,workPath) # Run time
            runTimeV[tIdx]= runTime
            tIdx = tIdx + 1
            print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                  np.mean(runTimeV), np.std(runTimeV)))
            print("--------------------------------------------------------")


    if(options == 'all' or options == 'cuffquant'):
        print("Quantifying gene expression using cuffquant")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        outDirPref = os.path.abspath("{}/output/rnaseq/".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        outDirPref = os.path.abspath("{}/output/rnaseq/cuffquant".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        inGtfDirPref  = os.path.abspath("{}/output/rnaseq/cuffmerge".format(workPath))   ## prefix
        inBamDirPref  = os.path.abspath("{}/output/rnaseq/tophat".format(workPath))   ## prefix
        ## Loop
        for size in rnaSeqSizeL:
            bamFileL  = glob.glob("{}/{}/*/accepted_hits.bam".format(inBamDirPref,size))
            outDir = "{}/{}".format(outDirPref,size)
            gtf="{}/{}/transcripts.gtf".format(inGtfDirPref,size)
            if(not os.path.isdir(outDir)):
                os.mkdir(outDir)

            for nThread in ompNumThreadsL:
                ## Consider adding nTrials here.
                runTimeV = np.zeros([len(bamFileL)])
                for bamFile in bamFileL:
                    outDirSamp = "{}/{}".format(outDir,bamFile.split("/")[-2].split(".")[0])
                    if(not os.path.isdir(outDirSamp)):
                        os.mkdir(outDirSamp)
                    tIdx = 0
                    cmd =  (
                        "time cuffquant --num-threads {} --output-dir {} "
                        "{} {}"
                            "".format(nThread, outDirSamp, gtf, bamFile))
                    output = "{}\n".format(cmd)
                    output = output + subprocess.getoutput(cmd)
                    runTime = parse_run_time(output,workPath) # Run time
                    runTimeV[tIdx]= runTime
                    tIdx = tIdx + 1
                print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
                print("--------------------------------------------------------")


    if(options == 'all' or options == 'cuffnorm'):
        print("Quantifying gene expression using cuffquant")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        outDirPref = os.path.abspath("{}/output/rnaseq/".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        outDirPref = os.path.abspath("{}/output/rnaseq/cuffnorm".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        inGtfDirPref  = os.path.abspath("{}/output/rnaseq/cuffmerge".format(workPath))   ## prefix
        inCxbDirPref  = os.path.abspath("{}/output/rnaseq/cuffquant".format(workPath))   ## prefix
        ## Loop
        for size in rnaSeqSizeL:
            cxbFileL  = glob.glob("{}/{}/*/abundances.cxb".format(inCxbDirPref,size))
            cxbFileL  = sorted(cxbFileL)    ## Break up into replicates
            # Get treat and wt groups
            sampNameL = [name.split('/')[-2] for name in cxbFileL]
            treatIdxL = ['treat_' in name for name in sampNameL]
            wtIdxL    = ['wt_' in name for name in sampNameL]
            treatCxbL = []
            wtCxbL    = []
            for idx in range(len(treatIdxL)):
                if(treatIdxL[idx] == True):
                    treatCxbL.append(cxbFileL[idx])
                elif(wtIdxL[idx] == True):
                    wtCxbL.append(cxbFileL[idx])
                else:
                    exit_with_error("ERROR!!! neither treatIdxL[idx] {} nor wtIdxL[idx] "
                                    "{} are" "True".format(treatIdxL[idx], wtIdxL[idx]))
            
            outDir = "{}/{}".format(outDirPref,size)
            gtf="{}/{}/transcripts.gtf".format(inGtfDirPref,size)
            if(not os.path.isdir(outDir)):
                os.mkdir(outDir)

            for nThread in ompNumThreadsL:
                ## Consider adding nTrials here.
                runTimeV = np.zeros([1])
                tIdx = 0
                cmd =  (
                    "time cuffnorm --num-threads {} --output-dir {} -L {} "
                      " {} {} {}"
                      "".format(nThread, outDir, "treat,wt",  gtf, 
                                ",".join(treatCxbL), ",".join(wtCxbL)))
                #print(cmd)
                output = "{}\n".format(cmd)
                output = output + subprocess.getoutput(cmd)
                runTime = parse_run_time(output,workPath) # Run time
                runTimeV[tIdx]= runTime
                tIdx = tIdx + 1
                print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
                print("--------------------------------------------------------")



    if(options == 'all' or options == 'cuffdiff'):
        print("Quantifying gene expression using cuffquant")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Size", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")
        outDirPref = os.path.abspath("{}/output/rnaseq/".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        outDirPref = os.path.abspath("{}/output/rnaseq/cuffdiff".format(workPath)) ## prefix
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        inGtfDirPref  = os.path.abspath("{}/output/rnaseq/cuffmerge".format(workPath))   ## prefix
        inCxbDirPref  = os.path.abspath("{}/output/rnaseq/cuffquant".format(workPath))   ## prefix
        ## Loop
        for size in rnaSeqSizeL:
            cxbFileL  = glob.glob("{}/{}/*/abundances.cxb".format(inCxbDirPref,size))
            cxbFileL  = sorted(cxbFileL)    ## Break up into replicates
            # Get treat and wt groups
            sampNameL = [name.split('/')[-2] for name in cxbFileL]
            treatIdxL = ['treat_' in name for name in sampNameL]
            wtIdxL    = ['wt_' in name for name in sampNameL]
            treatCxbL = []
            wtCxbL    = []
            for idx in range(len(treatIdxL)):
                if(treatIdxL[idx] == True):
                    treatCxbL.append(cxbFileL[idx])
                elif(wtIdxL[idx] == True):
                    wtCxbL.append(cxbFileL[idx])
                else:
                    exit_with_error("ERROR!!! neither treatIdxL[idx] {} nor wtIdxL[idx] "
                                    "{} are" "True".format(treatIdxL[idx], wtIdxL[idx]))
            
            outDir = "{}/{}".format(outDirPref,size)
            gtf="{}/{}/transcripts.gtf".format(inGtfDirPref,size)
            if(not os.path.isdir(outDir)):
                os.mkdir(outDir)

            for nThread in ompNumThreadsL:
                ## Consider adding nTrials here.
                runTimeV = np.zeros([1])
                tIdx = 0
                cmd =  (
                    "time cuffdiff --num-threads {} --output-dir {} -L {} "
                      " {} {} {}"
                      "".format(nThread, outDir, "treat,wt",  gtf, 
                                ",".join(treatCxbL), ",".join(wtCxbL)))
                #print(cmd)
                output = "{}\n".format(cmd)
                output = output + subprocess.getoutput(cmd)
                runTime = parse_run_time(output,workPath) # Run time
                runTimeV[tIdx]= runTime
                tIdx = tIdx + 1
                print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format(size, nThread,
                      np.mean(runTimeV), np.std(runTimeV)))
                print("--------------------------------------------------------")





    ## This will only run on Linux, not OSX
    if(options == 'all' or options == 'kelvin'):
        print("Runnning Kelving...")
        print("--------------------------------------------------------")
        print(" {:<10} | {:<12} | {:<15} | {:<15}".format("Short", "OMP_Threads", "mean",
              "stdev"))
        print("--------------------------------------------------------")

        # Create output directory structure
        outDirPref  = "{}/output".format(workPath)
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)
        outDirPref  = "{}/kelvin".format(outDirPref)
        curDir = os.path.dirname(os.path.realpath(__file__))
        if(not os.path.isdir(outDirPref)):
            os.mkdir(outDirPref)

        ## Loop
        for nThread in ompNumThreadsL:
            outDir = "{}/{}".format(outDirPref,nThread)
            if(not os.path.isdir(outDir)):
                os.mkdir(outDir)

            runTimeV = np.zeros([nTrials])
            for tIdx in range(nTrials):   ### change to nTrials
            ## Consider adding nTrials here.
                cmd =  (
                    "export OMP_NUM_THREADS={};"
                    "export LD_LIBRARY_PATH={}/src/kelvin/:$LD_LIBRARY_PATH;"
                    "time src/kelvin/kelvin src/kelvin/kelvin.conf --PedigreeFile src/kelvin/single.post > {}/kelvin.out.{}  2>&1"
                   "".format(nThread, curDir, outDir,tIdx))
                output = "{}\n".format(cmd)
                output = output + subprocess.getoutput(cmd)
                runTime = parse_run_time(output,workPath) # Run time
                runTimeV[tIdx]= runTime
            print(" {:<10} | {:<12} | {:<15.4f} | {:<15.4f}".format("Short", nThread,
                  np.mean(runTimeV), np.std(runTimeV)))
            print("--------------------------------------------------------")



    print("Run Time for {} option : {:.4f} h\n\n".format(options,(time.time() - startTime)/3600.0))
    sys.exit(0)


if __name__ == "__main__":
    main()
