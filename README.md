# Compute Node Benchmarking Tests
#### Author : Ali Snedden
#### Contact: Ali.Snedden[at]nationwidechildrens.org
#### License: MIT (unless noted otherwise, i.e. stream.c)
## Purpose:
We are seeking "Requests For Proposals" from various vendors of high performanch computing systems.
In order to make a selection from the vendors, we are providing small benchmarking test suite to be run on a typical compute node.
The tests are focused on utilities that commonly run on our compute current cluster. 
These include parts / programs used in bioinformatices, statistical genetics and mathematical modeling workflows.

The goal is to run each program multiple times, get a mean and standard deviation of the run time and to compare this for compute node various configurations from various vendors.
This Singularity recipe file should build everything you need within a CentOS 6.10 environment. 


## Installation:
#### Dependencies:
1. [Singularity](https://www.sylabs.io/guides/2.5/user-guide/index.html) (tested with Singularity 2.5.1)
2. A `/tmp` with about 20GB of free space.  
   This is used as a bind point for the Singularity image. 
   Most Linux machines should have this. This is used to place data and analysis output. 
3. A machine that you have root permission on. 

#### Build / Run:
Run : 
1.  Change the paths within the `%files` section of the `Singularity` recipe file to 
    the actual location of where the code lives.
    This is necessary because it is impossible to pass environmental variables to this section of the singularity recipe, which leads to this suboptimal solution:
    
    E.g.
    
    Change `/gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking/src` -> 
    `/your/unique/path/to/compute_node_benchmark/src`

    
2. `$ sudo singularity build test.simg Singularity`
3. `$ sudo chown user:group test.simg`
4. `$ singularity run -H /home/user test.simg`
5.  After inspecting the output files, `driver.log` (see `/tmp/benchmarking_out`) and stdout (from running the container), consider deleting `/tmp/benchmarking_out` to clean up your local disk.  
    `driver.log` contains the `stdout` from all the tests, but not `stdout` from running the container. 
    The run time values from the various programs will be printed to `stdout` when running step 3.
6. (If you are a vendor) Please return (via email) to us :
    * `driver.log`
    * The output from `stdout` of running the container (step 3).


#### Build Difficulties :
If you have difficulty with the downloads timing out or other network issues and the build (see above) is failing as a result, you can directly download the data into the cloned repo and then build the container.
E.g.

1. `$ cd /path/to/compute_node_benchmark`
2. `$ mkdir /ref`
3. `$ cd /opt/ref`
4. `$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz `
5. `$ tar xvzf grch38.tar.gz`
6. `$ mv grch38 HisatIndex`
7. `$ rm grch38.tar.gz`
8. `$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz`
9. `$ tar xvzf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz`
10. `$  mkdir Bowtie2Index`
11. `$  mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.1.bt2  Bowtie2Index/Homo_sapiens.GRC38.1.bt2`
12. `$  mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.2.bt2  Bowtie2Index/Homo_sapiens.GRC38.2.bt2`
13. `$  mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.3.bt2  Bowtie2Index/Homo_sapiens.GRC38.3.bt2`
14. `$  mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.4.bt2  Bowtie2Index/Homo_sapiens.GRC38.4.bt2`
15. `$  mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.1.bt2  Bowtie2Index/Homo_sapiens.GRC38.rev.1.bt2`
16. `$  mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.2.bt2  Bowtie2Index/Homo_sapiens.GRC38.rev.2.bt2`
17. `$  rm GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz`
18. `$  wget http://ftp.ensemblorg.ebi.ac.uk/pub/release-83/gtf/homo_sapiens/Homo_sapiens.GRCh38.83.gtf.gz`
19. `$  gunzip Homo_sapiens.GRCh38.83.gtf.gz`
20. `$  wget http://ftp.ensemblorg.ebi.ac.uk/pub/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
21. `$  gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
22. `$  cd ../`
23. Follow the steps in the Build / Run section

#### Warning:
1. When building the container it will be enormous (i.e. > 10GB). Beware!
2. The output is also enormous ~10-20GB. Beware!

#### Devs only:
1. If you want to unit test particular parts of the Singularity container, you can do this with something like :

   `$ export SINGULARITYENV_UNIT_TEST=cuffmerge`

   `$ export SINGULARITYENV_PREV_OUTPUT=/tmp/prev_benchmarking_out`

   `$ singularity run -H /home/group/user test.simg`

Note:
Typically one would copy `benchmarking_out` from a previous run to `/tmp/prev_benchmarking_out` and then `SINGULARITYENV_PREV_OUTPUT=/tmp/prev_benchmarking_out`


## What this code does:
1. Creates some matrix data (see `matrix_generator.py`)
2. Runs matrix multiplication tests with cache optimized code (see `matrix_multiply_omp_cache_optimized.c`). Tries various number of threads.
3. Runs matrix multiplication tests with non-cache optimized code (see `matrix_multiply_omp.c`)
4. Creates synthetic fastq data set (see `simulate_fastq.py`). Tries various number of threads.
5. Runs alignment to human genome using `tophat2`.
6. Runs alignment to human genome using `hisat2`.
7. Runs assembly on samples using `cufflinks`
8. Merges assembled data using `cuffmerge`
9. Runs `cuffcompare`
10. Quantifies gene/transcript expression with `cuffquant`
11. Normalizes gene/transcript expression with `cuffnorm`
12. Gets differential expression with `cuffdiff`
13. Runs `Kelvin` : which measures statistical evidence in human genetics
14. Runs `Stream` to test memory bandwidth


## Notes:
#### Development process:
This repository was extensively developed on OSX 10.13.6 using a python-3.6 virtual environment.
The `driver.py` (and other Python code that I developed) was run and tested primarily on this platform.
This saved time due how long it takes for the Singularity recipe to build.
After the code was in satisfactory shape, I then ported it over to run within a Singularity container.

The consequence of this development process is that there are some ugly portions of the code, where I'm:
1. Checking to see whether I'm using my local OSX python-3.6 virtual environment, or using the container's environment.
2. Within the recipe file I do a check of whether I have a local copy of the reference files (e.g. `bowtie`/`hisat`, genome and gtf files), if not I download directly. 
Downloading directly is slow, which is why I prefer the direct copy.
3. There are likely unforseen bugs due programs working on one platform (i.e. CentOS container or OSX).

#### Intended Machine to Test:
This is intended to test machine with : 
1. `x86_64` architecture, either AMD or Intel CPUs
2. By default we test up to 20 cores. We also test the number of cores equivalent to : a NUMA node, a socket and all the nodes cores.
3. These tests are not indended to test GPUs, FPGAs or the filesystem.

#### References:
Below are references for the code that we are testing. Most of this code (other than Stream and Kelvin) are downloaded from the internet at build time.
1. [Stream](http://www.cs.virginia.edu/stream/ref.html#why). See published [paper](https://www.researchgate.net/publication/51992086_Memory_bandwidth_and_machine_balance_in_high_performance_computers)
2. [Bowtie2 / Tophat2](https://ccb.jhu.edu/software/tophat/manual.shtml). See published [Paper](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-4-r36)
3. [HiSat2](https://ccb.jhu.edu/software/hisat2/index.shtml). See published [paper](https://www.nature.com/articles/nmeth.3317)
4. [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/index.html)
5. [Kelvin](https://www.karger.com/Article/Abstract/330634)



