BootStrap: yum
OSVersion: 6
MirrorURL: http://mirror.centos.org/centos-6/6/os/x86_64/
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
# Run : 
#   1. singularity run --nv -H /home/group/user test.simg   ## Runs prebuilt tests
#   2. singularity exec --nv -H /home/group/user test.simg  cmd  ## Runs cmd 
#


%help
Put Help here


#%environment

## Copy Files from outside of image to a place visible to image. Is Run before %post
%files
#/source /destination

#%labels
#    Author: Ali Snedden
#
# This is additional packages on top of the base image from MirrorURL
%post
    yum install python34
    yum install python34-libs
    yum install python34-numpy
    yum install python34-tools
    echo "Hello from inside the container"


%runscript
    ######## Matrix Multiply ########
    # Do compilation of files
    export OMP_NUM_THREADS=20
    gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp_cache_optimized.c -o src/matrix/matrix_multiply_omp_cache_optimized.o
    gcc -O3 -fopenmp src/matrix_multiply_omp.o -o src/matrix_multiply 

    # Generate files to run
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

    echo "Creating 10000 10000 10000 10000 files"
    mkdir -p data/10000/output
    python3 src/matrix/matrix_generator.py 10000 10000 10000 10000
    mv data/AB.txt data/10000
    mv data/B.txt  data/10000
    mv data/A.txt  data/10000

    
    echo "I've been run"
