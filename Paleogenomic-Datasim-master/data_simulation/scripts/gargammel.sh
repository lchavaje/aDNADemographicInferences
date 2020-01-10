#!/bin/bash

g_dir=$1
size_freq_path=$2
mat_file="/cm/shared/apps/gargammel/6oct2017/src/matrices/double-"

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

# Extract contamination and coverage values
coverages=( "${@:4:$3}" ); shift "$(( $3 + 3 ))"
contaminations=( "${@:2:$1}" ); shift "$(( $1 + 1 ))"

for cov in "${coverages[@]}"
do
	for cont in "${contaminations[@]}"
	do
		if [ "$cont" = "0" ]
		then
			complement=1
		else
			complement=0.$((100 - $cont))
		fi

		if [ "$cont" -gt "9" ]
		then
			per_cont="0."$cont
		else
			per_cont="0.0"$cont
		fi
		gargammel.pl -c "$cov" --comp 0,"$per_cont","$complement" -f "$size_freq_path" -matfile "$mat_file" -o "$g_dir"/cases/case_$SGE_TASK_ID/"$cov"x/"$cont"percent/frags "$g_dir"/gargammel/cases/case_$SGE_TASK_ID/ &
	done
done

# Copy original haplotypes to analyze later
cp "$g_dir"/gargammel/cases/case_$SGE_TASK_ID/endo/endo.1.fa "$g_dir"/cases/case_$SGE_TASK_ID/chr.1.fa
cp "$g_dir"/gargammel/cases/case_$SGE_TASK_ID/endo/endo.2.fa "$g_dir"/cases/case_$SGE_TASK_ID/chr.2.fa

wait
