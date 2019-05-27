# Compute Node Benchmarking Tests
#### Author : Ali Snedden
#### License: MIT (unless noted otherwise, i.e. stream.c)
## Purpose
We are seeking "Requests For Proposals" from various vendors of high performanch computing systems.
In order to make a selection from the vendors, we are providing small benchmarking test suite to be run on a typical compute node.
The tests are focused on utilities that commonly run on our compute current cluster. 
These include parts / programs used in bioinformatices, statistical genetics and mathematical modeling workflows.

The goal is to run each program multiple times, get a mean and standard deviation of the run time and to compare this for compute node various configurations from various vendors.
This Singularity recipe file should build everything you need within a CentOS 6.10 environment. 

## Notes
### Development process
This repository was extensively developed on OSX 10.13.6 using a python-3.6 virtual environment.
The `driver.py` (and other Python code that I developed) was run and tested primarily on this platform.
This saved time due how long it takes for the Singularity recipe to build.
After the code was in satisfactory shape, I then ported it over to run within a Singularity container.

The consequence of this development process is that there are some ugly portions of the code, where I'm:
1. Checking to see whether I'm using my local OSX python-3.6 virtual environment, or using the container's environment.
2. Within the recipe file I do a check of whether I have a local copy of the reference files (e.g. `bowtie`/`hisat`, genome and gtf files), if not I download directly. 
Downloading directly is slow, which is why I prefer the direct copy.
3. There are likely unforseen bugs due programs working on one platform (i.e. CentOS container or OSX).

### Intended Machine to Test
This is intended to test machine with : 
1. `x86_64` architecture, either AMD or Intel CPUs
2. Up to 20 cores are tested. More cores are acceptable, they just aren't tested in this code

These tests are not indended to test :
1. GPUs
2. FPGAs
3. The Filesystem.

### References
Below are references for the code that we are testing. Most of this code (other than Stream and Kelvin) are downloaded from the internet at build time.
1. [Stream](http://www.cs.virginia.edu/stream/ref.html#why). [Paper](https://www.researchgate.net/publication/51992086_Memory_bandwidth_and_machine_balance_in_high_performance_computers)
2. [Bowtie2 / Tophat2](https://ccb.jhu.edu/software/tophat/manual.shtml). [Paper](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-4-r36)
3. [HiSat2](https://ccb.jhu.edu/software/hisat2/index.shtml). [Paper](https://www.nature.com/articles/nmeth.3317)
4. [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/index.html)
5. [Kelvin](https://www.karger.com/Article/Abstract/330634)

## Installation
### Dependencies
1. Singularity (tested with Singularity 2.5.1)
2. A `/tmp` which is used as a bind point for the Singularity image. Most Linux machines should have this. This is used to place data and analysis output. 
3. A machine that you have root permission on. 

### Build / Run
Run : 
1. `sudo singularity build test.simg Singularity`
2. `sudo chown user:group test.simg`
3. `singularity run --nv -H /home/user test.simg`
4. After inspecting the output files, `driver.log` (see `/tmp/benchmarking_out`) and stdout (from running the container) , consider deleting `/tmp/benchmarking_out` to clean up your local disk.  `driver.log` contains the `stdout` from all the tests, but not `stdout` from running the container.  The run time values from the various programs will be printed to `stdout` when running step 3.

Note:
1. When building the container it will be enormous (i.e. > 10GB). Beware!

Please return to us :
1. `driver.log`
2. The output from `stdout` of running the container (step 3).
