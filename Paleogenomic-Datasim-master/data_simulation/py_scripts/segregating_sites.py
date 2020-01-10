#####################################################
# This script will find the segregating sites for   #
# all ancient human individuals in the './ancient'  #
# directory, gen_genomes.py should be run first     #
#####################################################
import sys
import os
from multiprocessing import Pool

# Global dir
g_dir = sys.argv[1]
anc_dir = g_dir + '/ancient/'

# Make sure directory exists
if not os.path.exists(anc_dir):
    print('segregating_sites.py: ancient sequence folder not present, '
            'try running gen_genomes.py first')
    quit()

# Get file list
files = [os.path.join(anc_dir, f) for f in os.listdir(anc_dir) if
        os.path.isfile(os.path.join(anc_dir, f))]

# Make sure it has the correct format, we need two chromosomes
# per individual
if len(files) % 2 != 0:
    print('segregating_sites.py: odd number of ancient chromosomes')
    quit()

# Calculate total individuals
individuals = int(len(files) / 2)

######################################
# Aux function, will be parallelized #
######################################

def write_sites(individual_index):
    # Individual string
    ind_string = 'ancient.' + str(individual_index + 1)

    # Build both filepaths
    chr_path_1 = anc_dir + ind_string + '.1.fa'
    chr_path_2 = anc_dir + ind_string + '.2.fa'

    # Segregating sites filepath
    seg_path = anc_dir + 'segsites.' + str(individual_index + 1)

    # Lines that will be written to segsites file
    out = []

    # Open both chromosome files
    with open(chr_path_1, 'r') as chr1, open(chr_path_2, 'r') as chr2:
        # Extract sequences
        line_1 = chr1.readlines()[1]
        line_2 = chr2.readlines()[1]

        # Compare base by base
        site = 1
        for x, y in zip(line_1, line_2):

            if x != y:
                out.append('>ref_1\t' + str(site) + '\t' + x.upper() +
                        '\t' + y.upper() + '\n')

            site += 1

    # Write to file
    with open(seg_path, 'w') as f:
        f.writelines(out)

#####################################################
# Iterate over chromosome files for all individuals #
#####################################################

# Parallelize over the number of individuals
pool = Pool()
pool.map(write_sites, range(individuals))
