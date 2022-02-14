import sys
import subprocess
import shutil
import os

########################
# Important parameters #
########################

# Global dir
g_dir = sys.argv[1]
# Number of cases, one for each ancient individual
num_cases = int(len(os.listdir(g_dir + '/present/')) / 3)
# Coverage parameters
#coverages = sys.argv[2].split(',')
# Contamination parameters
#contaminations = sys.argv[3].split(',')

#######################
# Directory structure #
#######################

# Directory constants
garg_dir = f'{g_dir}/gargammel/pre_cases/'
main_dir = f'{g_dir}/pre_cases/'

# Create clean directories, one for each case
if os.path.exists(main_dir):
    shutil.rmtree(main_dir)
if os.path.exists(garg_dir):
    shutil.rmtree(garg_dir)

os.makedirs(main_dir)
os.makedirs(garg_dir)

for i in range(num_cases):
    dir_string = main_dir + 'case_' + str(i + 1) + '/'
    garg_dir_string = garg_dir + 'case_' + str(i + 1) + '/'

    os.makedirs(dir_string)
    os.makedirs(garg_dir_string)

    # Directories for contaminant and endogenous data
    endo_dir = garg_dir_string + 'endo/'
    cont_dir = garg_dir_string + 'cont/'
    bact_dir = garg_dir_string + 'bact/'

    os.makedirs(endo_dir)
    os.makedirs(cont_dir)
    os.makedirs(bact_dir)

    # Directories for different coverage parameters
    #for c in coverages:
        #cov_dir = dir_string + c + 'x/'
        #os.makedirs(cov_dir)

        # Directories for different contamination parameters
        #for x in contaminations:
           #os.makedirs(cov_dir + x + 'percent/')

#################
# Data movement #
#################
#con_dir = g_dir + '/contaminant/'
ref_dir = g_dir + '/reference/'

# Make copies of contaminant and reference for each case
for i in range(num_cases):
    case_index = i + 1
    dir_string = garg_dir + 'case_' + str(case_index) + '/'

    #shutil.copyfile(con_dir + 'cont.1.fa', dir_string + 'cont/cont.1.fa')
    #shutil.copyfile(con_dir + 'cont.2.fa', dir_string + 'cont/cont.2.fa')

    shutil.copyfile(ref_dir + 'ref.fa', dir_string + 'ref.fa')

# Move each present human (and segsites) to its own case directory
for present in os.listdir(g_dir + '/present/'):
    filepath = g_dir + '/present/' + present
    case_index = present.split('.')[1]

    # Special case for segsites
    if 'segsites' in present:
        # Simply move it to its case folder
        dir_string = garg_dir + 'case_' + case_index + '/'
        shutil.move(filepath, dir_string + 'endo/segsites')
        continue

    # Get case number (individual's index) and chr index
    chr_index = present.split('.')[2]

    # Move file to case folder
    dir_string = garg_dir + 'case_' + case_index + '/'
    shutil.move(filepath, dir_string + 'endo/endo.' + str(chr_index) + '.fa')

###########
# Cleanup #
###########
if os.path.exists(g_dir + '/present'):
    shutil.rmtree(g_dir + '/present')

if os.path.exists(g_dir + '/contaminant'):
    shutil.rmtree(g_dir + '/contaminant')
