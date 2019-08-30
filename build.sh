#!/bin/bash
# Author : Ali Snedden
# Date : 8/29/19
# Purpose:
#   This script downloads and builds software / reference databases to be used
#   in run.sh
#
#   This will only work if you run this within the compute_node_benchmark directory
#   I use relative path names, so this is important
#
#
# HOW TO RUN:
#   sbatch --mail-type=FAIL,REQUEUE,TIME_LIMIT_90 --mail-user=XX@YY.com build.sh 
#
#
#
#
set -e


######## Reference Files ############
mkdir ref
cd ref
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

cd ../
cp src/simulate_fastq_data/data/chr1_short* ref

############# Do compilation of files ###################
## Cache optimized 
gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp_cache_optimized.c -o src/matrix/matrix_multiply_omp_cache_optimized.o
gcc -O3 -fopenmp src/matrix/matrix_multiply_omp_cache_optimized.o -o src/matrix/matrix_multiply_cache_opt
# Non-Cache optimized 
gcc -O3 -fopenmp -c src/matrix/matrix_multiply_omp.c -o src/matrix/matrix_multiply_omp.o
gcc -O3 -fopenmp src/matrix/matrix_multiply_omp.o -o src/matrix/matrix_multiply_non_cache_opt
## Compile stream with different array sizes
gcc -fopenmp -O -DSTREAM_ARRAY_SIZE=10000000 src/stream/stream.c -o src/stream/stream.10M
gcc -fopenmp -O -DSTREAM_ARRAY_SIZE=100000000 src/stream/stream.c -o src/stream/stream.100M


################ Install necessary software ################
mkdir software
cd    software
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
patch tophat-2.1.1 ../src/patch/tophat.patch

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
export PATH=`pwd`/software/samtools-1.9/bin:$PATH
cd ref
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
cd ../
echo "Finished building Software and Reference Data"


