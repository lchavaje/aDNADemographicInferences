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

coverages=$9
contaminations=${10}


module load hdf5/1.8.19
module load python37/3.7.0
module load seq-gen/1.3.4

echo "The variable 1 is equal to $1"
python3 py_scripts/Ms_parameters.py $1
i
# python3 $1/py_scripts/segregating_sites.py $1
# python3 $1/py_scripts/merge_reference.py $1
# python3 $1/py_scripts/gen_cases.py $1 $coverages $contaminations
# python3 $1/py_scripts/create_panel.py $1
