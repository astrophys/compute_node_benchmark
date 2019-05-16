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
#/source /destination

#%labels
#    Author: Ali Snedden
#
# This is additional packages on top of the base image from MirrorURL
%post
    yum install python34
    yum install python34-devel
    yum install python34-setuptools
    yum install python34-libs
    yum install python34-numpy
    yum install python34-tools

    mkdir -p /opt/python3.4/lib/python3.4/site-packages
    export PYTHONPATH=/opt/python3.4/lib/python3.4/site-packages:$PYTHONPATH
    export PYTHONPATH=/opt/python3.4/lib64/python3.4/site-packages:$PYTHONPATH
    export PATH=/opt/python3.4/bin:$PATH
    easy_install-3.4 --prefix /opt/python3.4 pip
    pip install --prefix /opt/python3.4 scipy
    
    #### Do compilation of files ####
    # Cache optimized 
    gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp_cache_optimized.c -o src/matrix/matrix_multiply_omp_cache_optimized.o
    gcc -O3 -fopenmp src/matrix/matrix_multiply_omp_cache_optimized.o -o src/matrix/matrix_multiply_cache_opt
    # Non-Cache optimized 
    gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp.c -o src/matrix/matrix_multiply_omp.o
    gcc -O3 -fopenmp src/matrix/matrix_multiply_omp.o -o src/matrix/matrix_multiply_non_cache_opt

    echo "Hello from inside the container"


%runscript
    ######## Matrix Multiply ########
    # Set Environmental Variables
    export PYTHONPATH=/opt/python3.4/lib/python3.4/site-packages:$PYTHONPATH

    # Generate files to run - This won't work here b/c it will try to write files to a directory and recall Singularity images cannot modify themselves in runscript. Recall we can use SINGULARITYENV_PATH
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

    #echo "Creating 10000 10000 10000 10000 files"
    #mkdir -p data/10000/output
    #python3 src/matrix/matrix_generator.py 10000 10000 10000 10000
    #mv data/AB.txt data/10000
    #mv data/B.txt  data/10000
    #mv data/A.txt  data/10000

    
    echo "I've been run"
