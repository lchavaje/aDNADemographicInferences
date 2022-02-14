#!/bin/bash

g_dir=$1
#size_freq_path=$2
#mat_file="/cm/shared/apps/gargammel/6oct2017/src/matrices/double-"

#$ -b y
#$ -w e
#$ -N gargammel_ref  # job's name
#$ -hold_jid pre_gargammel_ref

# run the job in the current working directory
#$ -cwd

# declares the shell that interprets the job script
#$ -S /bin/bash

# all environment variables in the qsub commands environment are to be exported to the batch job
#$ -V

# reserve all requested resources
#$ -l vf=8G
#$ -R y

module load htslib/1.2.1
module load samtools/1.2
module load seq-gen/1.3.4
module load gargammel/6oct2017


gargammel.pl -c 30 --comp 0,0,1 -l 150  -o "$g_dir"/pre_cases/case_$SGE_TASK_ID/frags "$g_dir"/gargammel/pre_cases/case_$SGE_TASK_ID/ &
	

# Copy original haplotypes to analyze later
cp "$g_dir"/gargammel/pre_cases/case_$SGE_TASK_ID/endo/endo.1.fa "$g_dir"/pre_cases/case_$SGE_TASK_ID/chr.1.fa
cp "$g_dir"/gargammel/pre_cases/case_$SGE_TASK_ID/endo/endo.2.fa "$g_dir"/pre_cases/case_$SGE_TASK_ID/chr.2.fa

wait
