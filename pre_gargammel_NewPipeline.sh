#!/bin/bash

#$ -b y
#$ -w e
#$ -N pre_gargammel_ref  # job's name

#$ -t 1
# run the job in the current working directory
#$ -cwd

# declares the shell that interprets the job script
#$ -S /bin/bash

# all environment variables in the qsub commands environment are to be exported to the batch job
#$ -V

# reserve all requested resources
#$ -l vf=70G
#$ -R y

coverages=$3
contaminations=$4

echo "coverages is equal to $3 and contaminations is equal to $4"

module load hdf5/1.8.19
module load python37/3.7.0
module load seq-gen/1.3.4

#echo "The variable 1 is equal to $1 and the variable 2 is equal to $2"
python3 py_scripts/sample_simulation.py $1 $2
python3 $2/py_scripts/anc_seg_sites.py $2
echo "anc_seg_sites is finished"
python3 $2/py_scripts/merge_reference.py $2
echo "merge_reference is finished"
python3 $2/py_scripts/pre_seg_sites.py $2
echo "pre_seg_sites is finished"
python3 $2/py_scripts/gen_anc_cases.py $2 $coverages $contaminations
echo "gen_anc_cases is finished"
python3 $2/py_scripts/gen_pre_cases.py $2 $coverages $contaminations
echo "gen_pre_cases is finished"
#python3 $2/py_scripts/create_panel.py $2
#echo "create_panel is finished"
