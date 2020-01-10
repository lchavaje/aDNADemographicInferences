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
module load python36/3.6.3
module load seq-gen/1.3.4

if [ "$7" = "null" ] || [ "$8" = "null" ]
then
	python3 $1/py_scripts/gen_genomes.py -d $1 -g $2 -a $3 -r $4 -b $5 -s $6
else
	python3 $1/py_scripts/gen_genomes.py -d $1 -g $2 -a $3 -r $4 -b $5 -s $6 -e $7 -t $8
fi

python3 $1/py_scripts/segregating_sites.py $1
python3 $1/py_scripts/merge_reference.py $1
python3 $1/py_scripts/gen_cases.py $1 $coverages $contaminations
python3 $1/py_scripts/create_panel.py $1
