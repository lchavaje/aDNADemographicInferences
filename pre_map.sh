#!/bin/bash

# Declare bash variables
g_dir=$1
ref=$1"/reference/ref.fa"

#$ -b y
#$ -w e
#$ -N pre_map_ref  # job's name
#$ -hold_jid gargammel_ref

#$ -t 1
# run the job in the current working directory
#$ -cwd

# declares the shell that interprets the job script
#$ -S /bin/bash

# all environment variables in the qsub commands environment are to be exported to the batch job
#$ -V

# reserve all requested resources
#$ -l vf=8G
#$ -R y

# Cleanup
rm -rf "$g_dir"/gargammel/

# Index reference genome
module load htslib/1.2.1
module load bwa/0.7.15
module load samtools/1.2

bwa index -a bwtsw "$ref"
samtools faidx "$ref"
