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

    echo "Hello from inside the container"


%runscript
    ######## Matrix Multiply ########
    # Do compilation of files
    export OMP_NUM_THREADS=20
    gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp_cache_optimized.c -o src/matrix/matrix_multiply_omp_cache_optimized.o
    gcc src/matrix_multiply_omp.o -o src/matrix_multiply 

    # Generate files to run
    
    echo "I've been run"
