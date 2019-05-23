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
#      sudo SINGULARITYENV_CODEPATH=/some/path/to/the/benchmarking_dir; singularity build test.simg Singularity
#      sudo chown user:group test.simg
# Run : 
#   0. Turn off multithreading
#   1. singularity run --nv -H /home/group/user test.simg   ## Runs prebuilt tests
#   2. singularity exec --nv -H /home/group/user test.simg  cmd  ## Runs cmd 
#


%help
Put Help here


#%environment

## Copy Files from outside of image to a place visible to image. Is Run before %post
%files
# Will need to use environmental variables to copy the code to 
/gpfs0/home/gdhpcgroup/aps003/Code/Singularity/benchmarking /tmp
/reference/homo_sapiens/GRCh38/ensembl/release-83/Annotation/Genes/gtf/Homo_sapiens.GRCh38.83.gtf /tmp
/reference/homo_sapiens/GRCh38/ensembl/release-83/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa /tmp

#%labels
#    Author: Ali Snedden
#
# This is additional packages on top of the base image from MirrorURL
%post
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

    mkdir -p /opt/python3.4/lib/python3.4/site-packages
    export PYTHONPATH=/opt/python3.4/lib/python3.4/site-packages:$PYTHONPATH
    export PYTHONPATH=/opt/python3.4/lib64/python3.4/site-packages:$PYTHONPATH
    export PATH=/opt/python3.4/bin:$PATH
    easy_install-3.4 --prefix /opt/python3.4 pip
    pip install --prefix /opt/python3.4 scipy
    #
    ##### Do compilation of files ####
    mkdir /opt/code
    mv /tmp/benchmarking /opt/code/

    ## Cache optimized 
    cd /opt/code/benchmarking
    gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp_cache_optimized.c -o src/matrix/matrix_multiply_omp_cache_optimized.o
    gcc -O3 -fopenmp src/matrix/matrix_multiply_omp_cache_optimized.o -o src/matrix/matrix_multiply_cache_opt
    # Non-Cache optimized 
    gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp.c -o src/matrix/matrix_multiply_omp.o
    gcc -O3 -fopenmp src/matrix/matrix_multiply_omp.o -o src/matrix/matrix_multiply_non_cache_opt

    # Get genomics data
    mkdir /opt/ref
    mv /tmp/*.fa  /opt/ref
    mv /tmp/*.gtf /opt/ref
    ### Uncomment if this isn't locally located ###
    #wget http://ftp.ensemblorg.ebi.ac.uk/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz
    #gunzip Homo_sapiens.GRCh38.96.gtf.gz
    #wget http://ftp.ensemblorg.ebi.ac.uk/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    #gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

    # Install necessary software
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

    echo "Hello from inside the container"


%runscript
    ######## Matrix Multiply ########
    # Set Environmental Variables
    export PYTHONPATH=/opt/python3.4/lib/python3.4/site-packages:$PYTHONPATH
    export PATH=/opt/software/tophat-2.1.1:$PATH
    export PATH=/opt/software/hisat2-2.1.0:$PATH
    export PATH=/opt/software/bowtie2-2.3.5.1:$PATH

    # Generate files to run - This won't work here b/c it will try to write files to a directory and recall Singularity images cannot modify themselves in runscript. Recall we can use SINGULARITYENV_PATH
    ## I think python3 src/driver.py all will replace below code - be sure it writes, reads, and cleans up tmp 
    echo "Creating 2000 10000 10000 2000 files"
    mkdir -p data/2000/output
    python3 src/matrix/matrix_generator.py 2000 10000 10000 2000
    mv data/AB.txt data/2000
    mv data/B.txt  data/2000
    mv data/A.txt  data/2000

    echo "Creating 3000 10000 10000 3000 files"
    mkdir -p data/3000/output
    python3 src/matrix/matrix_generator.py 3000 10000 10000 3000
    mv data/AB.txt data/3000
    mv data/B.txt  data/3000
    mv data/A.txt  data/3000

    echo "Creating 4000 10000 10000 4000 files"
    mkdir -p data/4000/output
    python3 src/matrix/matrix_generator.py 4000 10000 10000 4000
    mv data/AB.txt data/4000
    mv data/B.txt  data/4000
    mv data/A.txt  data/4000

    echo "Creating 5000 10000 10000 5000 files"
    mkdir -p data/5000/output
    python3 src/matrix/matrix_generator.py 5000 10000 10000 5000
    mv data/AB.txt data/5000
    mv data/B.txt  data/5000
    mv data/A.txt  data/5000

    mkdir -p data/rnaseq/fastq
    #echo "Creating 10000 10000 10000 10000 files"
    #mkdir -p data/10000/output
    #python3 src/matrix/matrix_generator.py 10000 10000 10000 10000
    #mv data/AB.txt data/10000
    #mv data/B.txt  data/10000
    #mv data/A.txt  data/10000
    

    ### Run stream multiple times ###

    
    echo "I've been run"
