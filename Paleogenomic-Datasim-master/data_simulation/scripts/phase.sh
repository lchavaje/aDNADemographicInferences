#!/bin/bash

# Declare bash variables
g_dir=$1
ref_panel_hap="$g_dir""/reference/ref_panel.hap"
ref_panel_leg="$g_dir""/reference/ref_panel.leg"
ref_panel_sam="$g_dir""/reference/ref_panel.sam"

#$ -b y
#$ -w e
#$ -N phase_ref  # job's name
#$ -hold_jid call_variants_ref

# run the job in the current working directory
#$ -cwd

# declares the shell that interprets the job script
#$ -S /bin/bash

# all environment variables in the qsub commands environment are to be exported to the batch job
#$ -V

# reserve all requested resources
#$ -l vf=16G
#$ -R y

module load shapeit/2r837
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
	mkdir phased
	mkdir phased/logs

	#########
	# Phase #
	#########

	# Get sites to exclude because of misalignment
	shapeit -check \
		--input-vcf variants_biallelic.vcf.gz \
		--input-ref "$ref_panel_hap" "$ref_panel_leg" "$ref_panel_sam" \
		--output-log readset.alignments

	# Phase, ignoring sites in exclude file
	if [ -f readset.alignments.snp.strand.exclude ]; then
		shapeit --input-vcf variants_biallelic.vcf.gz \
			--input-ref "$ref_panel_hap" "$ref_panel_leg" "$ref_panel_sam" \
			--exclude-snp readset.alignments.snp.strand.exclude \
			--rho 0.00000002 \
			-O readset.phased	
	else
		shapeit --input-vcf variants_biallelic.vcf.gz \
			--input-ref "$ref_panel_hap" "$ref_panel_leg" "$ref_panel_sam" \
			--rho 0.00000002 \
			-O readset.phased	
	fi

	# Cleanup
	mv readset.alignments* phased/
	mv readset.phased* phased/
	mv shapeit* phased/logs

	# Get results
	python3 "$g_dir"/py_scripts/get_accuracy.py
	python3 "$g_dir"/py_scripts/get_switch_err.py
done
