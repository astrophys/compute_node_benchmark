BootStrap: yum
OSVersion: 6
MirrorURL: http://mirror.centos.org/centos-6/6/os/x86_64/
Include: yum
# If you want the updates (available at the bootstrap date) to be installed
# inside the container during the bootstrap instead of the General Availability
# point release (7.x) then uncomment the following line
# UpdateURL: http://mirror.centos.org/centos-%{OSVERSION}/%{OSVERSION}/updates/$basearch/

# Notes : 
#
# Build : 
#   1. rm test.simg
#      sudo singularity build test.simg Singularity
#      sudo chown user:group test.simg
# Run (all tests) : 
#   0. Turn off multithreading
#   1. singularity run --nv -H /home/group/user test.simg   ## Runs prebuilt tests
#
# Run (unit_test) : 
#   0. Turn off multithreading
#   
#   1. E.g. to test cuffmerge option to driver.py :
#            export SINGULARITYENV_UNIT_TEST=cuffmerge\
#            export SINGULARITYENV_PREV_OUTPUT=/tmp/prev_benchmarking_out/\
#            singularity run --nv -H /home/group/user test.simg   
#      --> Typically one would copy benchmarking_out from a previous run to 
#          /tmp/benchmarking_out and then 
#      --> Set SINGULARITYENV_PREV_OUTPUT=/tmp/benchmarking_out
#
#


%help
Put Help here


#%environment

## Copy Files from outside of image to a place visible to image. Is Run before %post
%files
# Will need to use environmental variables to copy the code to 
#/gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking/ /tmp  # Copies *.simg, maybe comment below later and uncomment this line?
/gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking/src /tmp
/gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking/WHERE_I_LEFT_OFF.txt /tmp
/gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking/config /tmp
/gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking/README.md /tmp
/gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking/Singularity tmp
/gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking/ref tmp

#%labels
#    Author: Ali Snedden
#
# This is additional packages on top of the base image from MirrorURL
%post
    echo "CODEPATH=$CODEPATH"
    ls /usr/bin
    yum install -y yum-utils.noarch
    yum install -y epel-release.noarch
    yum install -y python34
    yum install -y python34-devel
    yum install -y python34-setuptools
    yum install -y python34-libs
    yum install -y python34-numpy
    yum install -y python34-tools
    yum install -y gcc
    yum install -y wget
    yum install -y tar.x86_64
    yum install -y unzip.x86_64
    yum install -y autoconf.noarch
    yum install -y automake.noarch
    yum install -y time.x86_64
    yum install -y time.x86_64
    yum install -y util-linux-ng-2.17.2-12.28.el6_9.2.x86_64
    # Install python 2
    yum install -y python
    yum install -y python-devel
    yum install -y python-setuptools
    yum install -y python-libs
    #yum install -y python-numpy
    yum install -y python-tools
    
    # Samtools dependencies here
    yum install -y bzip2
    yum install -y bzip2-devel
    yum install -y zlib.x86_64
    yum install -y zlib-devel.x86_64
    yum install -y gawk.x86_64
    yum install -y lzip.x86_64
    yum install -y xz-devel.x86_64

    # Kelvin dependencies
    yum install -y gsl.x86_64
    yum install -y gsl-devel.x86_64
    

    mkdir -p /opt/python3.4/lib/python3.4/site-packages
    export PYTHONPATH=/opt/python3.4/lib/python3.4/site-packages:$PYTHONPATH
    export PYTHONPATH=/opt/python3.4/lib64/python3.4/site-packages:$PYTHONPATH
    export PATH=/opt/python3.4/bin:$PATH
    easy_install-3.4 --prefix /opt/python3.4 pip
    pip install --prefix /opt/python3.4 scipy
    pip install --prefix /opt/python3.4 argparse
    #
    ##### Do compilation of files ####
    mkdir -p /opt/code/benchmarking
    mv /tmp/src /opt/code/benchmarking
    mv /tmp/WHERE_I_LEFT_OFF.txt /opt/code/benchmarking
    mv /tmp/config /opt/code/benchmarking
    mv /tmp/README.md /opt/code/benchmarking
    mv /tmp/Singularity /opt/code/benchmarking
    mv /tmp/ref /opt/code/benchmarking

    # If reference directory exists, use it
    if [ -d "/opt/code/benchmarking/ref" ]; then 
        mv /opt/code/benchmarking/ref /opt
    # If not, download
    else
        ## Untested
        mkdir /opt/ref
        cd /opt/ref
        # Hisat indices
        wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
        tar xvzf grch38.tar.gz
        mv grch38 HisatIndex
        rm grch38.tar.gz
        # Bowtie indices
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
        tar xvzf GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
        mkdir Bowtie2Index
        mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.1.bt2  Bowtie2Index/Homo_sapiens.GRC38.1.bt2
        mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.2.bt2  Bowtie2Index/Homo_sapiens.GRC38.2.bt2
        mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.3.bt2  Bowtie2Index/Homo_sapiens.GRC38.3.bt2
        mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.4.bt2  Bowtie2Index/Homo_sapiens.GRC38.4.bt2
        mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.1.bt2  Bowtie2Index/Homo_sapiens.GRC38.rev.1.bt2
        mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.2.bt2  Bowtie2Index/Homo_sapiens.GRC38.rev.2.bt2
        rm GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz
        # Gtf
        wget http://ftp.ensemblorg.ebi.ac.uk/pub/release-83/gtf/homo_sapiens/Homo_sapiens.GRCh38.83.gtf.gz
        gunzip Homo_sapiens.GRCh38.83.gtf.gz
        # Fasta
        wget http://ftp.ensemblorg.ebi.ac.uk/pub/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    fi
    mv /opt/code/benchmarking/src/simulate_fastq_data/data/chr1_short* /opt/ref

    ##### Do compilation of files ####
    ## Cache optimized 
    cd /opt/code/benchmarking
    gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp_cache_optimized.c -o src/matrix/matrix_multiply_omp_cache_optimized.o
    gcc -O3 -fopenmp src/matrix/matrix_multiply_omp_cache_optimized.o -o src/matrix/matrix_multiply_cache_opt
    # Non-Cache optimized 
    gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp.c -o src/matrix/matrix_multiply_omp.o
    gcc -O3 -fopenmp src/matrix/matrix_multiply_omp.o -o src/matrix/matrix_multiply_non_cache_opt
    ## Compile stream with different array sizes
    gcc -fopenmp -O -DSTREAM_ARRAY_SIZE=10000000 src/stream/stream.c -o src/stream/stream.10M
    gcc -fopenmp -O -DSTREAM_ARRAY_SIZE=100000000 src/stream/stream.c -o src/stream/stream.100M

    ###### Install necessary software ######
    mkdir /opt/software
    cd    /opt/software
    # Bowtie
    wget https://github.com/BenLangmead/bowtie2/releases/download/v2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip
    unzip bowtie2-2.3.5.1-linux-x86_64.zip
    rm bowtie2-2.3.5.1-linux-x86_64.zip
    mv bowtie2-2.3.5.1-linux-x86_64 bowtie2-2.3.5.1

    # Tophat
    wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
    tar xvzf tophat-2.1.1.Linux_x86_64.tar.gz
    rm tophat-2.1.1.Linux_x86_64.tar.gz
    mv tophat-2.1.1.Linux_x86_64 tophat-2.1.1

    # Hisat
    wget http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-Linux_x86_64.zip
    unzip hisat2-2.1.0-Linux_x86_64.zip
    rm hisat2-2.1.0-Linux_x86_64.zip

    # Cufflinks
    wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
    tar xvzf cufflinks-2.2.1.Linux_x86_64.tar.gz
    rm cufflinks-2.2.1.Linux_x86_64.tar.gz
    mv cufflinks-2.2.1.Linux_x86_64 cufflinks-2.2.1
    
    # Samtools
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar xvjf samtools-1.9.tar.bz2
    rm samtools-1.9.tar.bz2
    cd samtools-1.9
    ./configure --without-curses
    make
    mkdir bin
    cp samtools bin/
    # Create fai for cufflinks workflow
    export PATH=/opt/software/samtools-1.9/bin:$PATH
    cd /opt/ref
    if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai" ]; then 
        samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
    fi

    # String Tie
    wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.6.Linux_x86_64.tar.gz
    tar xvzf stringtie-1.3.6.Linux_x86_64.tar.gz
    rm stringtie-1.3.6.Linux_x86_64.tar.gz
    mv stringtie-1.3.6.Linux_x86_64 stringtie-1.3.6

    # Ballgown
    # Must download via bioconductor : 
    # source("http://bioconductor.org/biocLite.R")
    # biocLite("ballgown")


    ########### Get genomics data ##########

    echo "Hello from inside the container"


%runscript
    set -e
    ######## Matrix Multiply ########
    # Set Environmental Variables
    export PYTHONPATH=/opt/python3.4/lib/python3.4/site-packages:$PYTHONPATH
    export PATH=/opt/software/tophat-2.1.1:$PATH
    export PATH=/opt/software/hisat2-2.1.0:$PATH
    export PATH=/opt/software/bowtie2-2.3.5.1:$PATH
    export PATH=/opt/software/cufflinks-2.2.1:$PATH
    export PATH=/opt/software/samtools-1.9/bin:$PATH
    export PATH=/opt/software/stringtie-1.3.6:$PATH

    ## Unit Tests ##
    if [[ ! -z "${UNIT_TEST}" ]]; then
        cd /opt/code/benchmarking
        echo "Running Unit Test : ${UNIT_TEST}"
        ## For most of these unit tests, previous output must have been generated
        ## so we'll dump the current output to the same dir as PREV_OUTPUT
        if [[ -z "${PREV_OUTPUT}" ]]; then
            echo "ERROR!!! UNIT_TEST set but PREV_OUTPUT is unset" >&2
            exit 1
        fi

        rm ${PREV_OUTPUT}/driver.log

        ## Error check that a correct option 
        if [ ${UNIT_TEST} == "build_mat_mult_data" ]; then 
            python3 src/driver.py build_mat_mult_data ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "mat_mult_cache_opt" ]; then 
            python3 src/driver.py mat_mult_cache_opt ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "mat_mult_non_cache_opt" ]; then 
            python3 src/driver.py mat_mult_non_cache_opt ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "build_rnaseq_data" ]; then 
            python3 src/driver.py build_rnaseq_data ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "align_rnaseq_tophat" ]; then 
            python3 src/driver.py align_rnaseq_tophat ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "align_rnaseq_hisat" ]; then 
            python3 src/driver.py align_rnaseq_hisat ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "cufflinks_assemble" ]; then 
            python3 src/driver.py cufflinks_assemble ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "cuffmerge" ]; then 
            python3 src/driver.py cuffmerge ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "cuffcompare" ]; then 
            python3 src/driver.py cuffcompare ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "cuffquant" ]; then 
            python3 src/driver.py cuffquant ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "cuffnorm" ]; then 
            python3 src/driver.py cuffnorm ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "cuffdiff" ]; then 
            python3 src/driver.py cuffdiff ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "kelvin" ]; then 
            python3 src/driver.py kelvin ${PREV_OUTPUT} /opt/ref/
        elif [ ${UNIT_TEST} == "stream" ]; then 
            echo "Running STREAM_ARRAY_SIZE = 10M" 
            echo "Run 1"
            bash src/stream/run_stream.sh src/stream/stream.10M
            echo "Run 2"
            bash src/stream/run_stream.sh src/stream/stream.10M
            echo "Run 3"
            bash src/stream/run_stream.sh src/stream/stream.10M
            echo "Run 4"
            bash src/stream/run_stream.sh src/stream/stream.10M
            echo "Run 5"
            bash src/stream/run_stream.sh src/stream/stream.10M
            echo "Running STREAM_ARRAY_SIZE = 100M" 
            echo "Run 1"
            bash src/stream/run_stream.sh src/stream/stream.100M
            echo "Run 2"
            bash src/stream/run_stream.sh src/stream/stream.100M
            echo "Run 3"
            bash src/stream/run_stream.sh src/stream/stream.100M
            echo "Run 4"
            bash src/stream/run_stream.sh src/stream/stream.100M
            echo "Run 5"
            bash src/stream/run_stream.sh src/stream/stream.100M
        else
            echo "ERROR!!! UNIT_TEST = ${UNIT_TEST} is an invalid option"
            exit 1
        fi

    ## Full test ##
    else

        if [ -d "/tmp/benchmarking_out" ]; then 
            rm -r /tmp/benchmarking_out
        fi
        mkdir -p /tmp/benchmarking_out
        mkdir -p /tmp/benchmarking_out/data
        mkdir -p /tmp/benchmarking_out/output
        cd /opt/code/benchmarking
        ## Temporarly comment out
        python3 src/driver.py build_mat_mult_data /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py mat_mult_cache_opt /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py mat_mult_non_cache_opt /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py build_rnaseq_data /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py align_rnaseq_tophat /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py align_rnaseq_hisat /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py cufflinks_assemble /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py cuffmerge /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py cuffcompare /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py cuffquant /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py cuffnorm /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py cuffdiff /tmp/benchmarking_out/ /opt/ref/
        python3 src/driver.py kelvin /tmp/benchmarking_out/ /opt/ref/


        ### Run stream multiple times ###
        echo "Running STREAM_ARRAY_SIZE = 10M" 
        echo "Run 1"
        bash src/stream/run_stream.sh src/stream/stream.10M
        echo "Run 2"
        bash src/stream/run_stream.sh src/stream/stream.10M
        echo "Run 3"
        bash src/stream/run_stream.sh src/stream/stream.10M
        echo "Run 4"
        bash src/stream/run_stream.sh src/stream/stream.10M
        echo "Run 5"
        bash src/stream/run_stream.sh src/stream/stream.10M
        echo "Running STREAM_ARRAY_SIZE = 100M" 
        echo "Run 1"
        bash src/stream/run_stream.sh src/stream/stream.100M
        echo "Run 2"
        bash src/stream/run_stream.sh src/stream/stream.100M
        echo "Run 3"
        bash src/stream/run_stream.sh src/stream/stream.100M
        echo "Run 4"
        bash src/stream/run_stream.sh src/stream/stream.100M
        echo "Run 5"
        bash src/stream/run_stream.sh src/stream/stream.100M
    fi

    
    echo "WARNING!!!" 
    echo "  After inspecting the containts of /tmp/benchmarking_out, you may "
    echo "  consider deleteing it"
