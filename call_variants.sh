#!/bin/bash

# Declare bash variables
g_dir=$1
ref=$g_dir"/reference/ref.fa"

#$ -b y
#$ -w e
#$ -N call_variants_ref  # job's name
#$ -hold_jid map_ref

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
module load samtools/1.6
module load bcftools/1.2
module load vcftools/0.1.14

module load hdf5/1.8.19
module load python36/3.6.3

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

	#################
	# Call variants #
	#################
	samtools mpileup -t AD,INFO/AD -B -v -f "$ref" readset_sort_nodup.bam | bcftools call -v -m -f GQ -O z - > variants.vcf.gz
	vcftools --gzvcf variants.vcf.gz --minQ 90 --remove-indels --recode --stdout | gzip - > variants_filtered.vcf.gz
	vcftools --gzvcf variants_filtered.vcf.gz --min-alleles 1 --max-alleles 2 --recode --stdout | gzip - > variants_biallelic.vcf.gz

	rm variants.vcf.gz
	rm variants_filtered.vcf.gz

	#python3 "$g_dir"/py_scripts/real_vcf.py $g_dir
	zcat variants_biallelic.vcf.gz > variants_biallelic.vcf
	python3 "$g_dir"/py_scripts/vcf_stats.py $g_dir
	rm variants_biallelic.vcf

	##########################
	# Store pre-filter depth #
	##########################
	samtools depth -a readset_sort_nodup.bam | awk '{sum+=$3} END {print sum/NR}' > depth
	samtools depth -a readset_sort_nodup.bam | awk '$3==0 {print $2}' > not_covered
done
