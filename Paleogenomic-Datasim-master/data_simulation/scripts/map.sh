#!/bin/bash

# Declare bash variables
g_dir=$1
ref=$1"/reference/ref.fa"

#$ -b y
#$ -w e
#$ -N map_ref  # job's name
#$ -hold_jid pre_map_ref

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
module load fastqc/0.11.3
module load trimgalore/0.4.2
module load bwa/0.7.15
module load samtools/1.2

# Extract contamination and coverage values
coverages=( "${@:3:$2}" ); shift "$(( $2 + 2 ))"
contaminations=( "${@:2:$1}" ); shift "$(( $1 + 1 ))"

# Transform into dirs
dirs=()
for cov in "${coverages[@]}"
do
	for cont in "${contaminations[@]}"
	do
		dirs+=("$g_dir""/cases/case_"$SGE_TASK_ID"/"$cov"x/"$cont"percent/")
	done
done

for dir in "${dirs[@]}"
do
	cd "$dir"

	###################
	# Adapter removal #
	###################

	# Extract compressed reads
	gunzip frags_s1.fq.gz
	gunzip frags_s2.fq.gz

	# Trim reads
	trim_galore --paired --length 25 frags_s1.fq frags_s2.fq

	# Cleanup
	mv *val_1.fq clean_R1.fastq
	mv *val_2.fq clean_R2.fastq
	rm frags*

	###########
	# Mapping #
	###########

	# BWA assembly
	bwa mem "$ref" clean_R1.fastq clean_R2.fastq > readset.sam
	samtools view -S readset.sam -b -o readset.bam

	# Sort and remove duplicates
	samtools sort readset.bam readset_sort
	samtools rmdup -s readset_sort.bam readset_sort_nodup.bam

	# Create index
	samtools index readset_sort_nodup.bam

	# Cleanup
	rm clean_R*
	rm readset.bam readset.sam readset_sort.bam
done
