#####################################################
# This script will merge two chromosomes to use as  #
# a mapping reference                               #
#####################################################
import sys
import os
import random

# Global dir
g_dir = sys.argv[1]
pre_dir = f'{g_dir}/present/'

# Make sure directory exists
if not os.path.exists(pre_dir):
    print('merge_reference.py: present sequence folder not present, '
            'try running ample_simulation.py first')
    quit()

# Get file list
files = [os.path.join(pre_dir, f) for f in os.listdir(pre_dir) if
        os.path.isfile(os.path.join(pre_dir, f))]

# Make sure it has the correct format, we need two chromosomes
# per individual
if len(files) % 2 != 0 or len(files) < 2:
    print('merge_reference.py: odd number of present chromosomes or not enough chromosomes')
    quit()

# Calculate total individuals
individuals = int(len(files) / 2)

# Individual string
ind_string = 'present.' + str(individuals)

# Build both filepaths
chr_path_1 = pre_dir + ind_string + '.1.fa'
chr_path_2 = pre_dir + ind_string + '.2.fa'

# Segregating sites filepath
seg_path = pre_dir + 'ref.fa'

# Lines that will be written to segsites file
out = []
out.append('>ref_1\n')

# Open both chromosome files
with open(chr_path_1, 'r') as chr1, open(chr_path_2, 'r') as chr2:
    # Extract sequences
    line_1 = chr1.readlines()[1]
    line_2 = chr2.readlines()[1]

    # Compare base by base
    for x, y in zip(line_1, line_2):

        if x != y:
                if random.random() < 0.5:
                    out.append(x.upper())
                else:
                    out.append(y.upper())
        else:
            out.append(x.upper())

# Write to file
with open(seg_path, 'w') as f:
    f.writelines(out)

os.rename(chr_path_1, g_dir + '/reference/chr1.fa')
os.rename(chr_path_2, g_dir + '/reference/chr2.fa')
os.rename(seg_path, g_dir + '/reference/ref.fa')
