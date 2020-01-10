#!/bin/bash

g_dir=$1

#$ -b y
#$ -w e
#$ -N results_ref  # job's name
#$ -hold_jid phase_ref

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

module load hdf5/1.8.19
module load python36/3.6.3

python3 "$g_dir"/py_scripts/get_results.py $1 $2 $3 $4 $5

mkdir "$g_dir"/"$2"gen/
mv "$g_dir"/cases "$g_dir"/"$2"gen/
mv "$g_dir"/reference "$g_dir"/"$2"gen/
mv "$g_dir"/present "$g_dir"/"$2"gen/
mv "$g_dir"/gargammel* "$g_dir"/"$2"gen/
