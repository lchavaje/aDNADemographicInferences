import os
import sys
import random

######################################################
# This script will find create a reference panel     #
# based on the individuals in the present/ directory #
######################################################

# Global dir
g_dir = sys.argv[1]
out_dir = f'{g_dir}/reference/'
pre_dir = f'{g_dir}/present/'
segsite_filepath = f'{g_dir}/all_segsites'

# Lists where file output will be stored, included headers
sam_f = ['sample population group sex']
leg_f = ['id position a0 a1']
hap_f = []

# Make sure directory exists
if not os.path.exists(pre_dir):
    print('create_panel.py: present sequence folder not present, '
            'try running gen_genomes.py first')
    quit()

if not os.path.isfile(segsite_filepath):
    print('create_panel.py: segsite file missing')
    quit()

# Get file list
files = [os.path.join(pre_dir, f) for f in os.listdir(pre_dir) if
        os.path.isfile(os.path.join(pre_dir, f))]

# Make sure it has the correct format, we need two chromosomes
# per individual, plus segsites
if len(files) % 2 != 0:
    print('create_panel.py: odd number of present chromosomes')
    quit()

# Calculate total individuals
individuals = int(len(files) / 2)

# List where all haplotypes will be stored
haplotypes = []

for i in range(individuals):
    # Add to SAMPLE file, dummy values for population
    # and group, random sex
    sam_f.append('IND' + str(i + 1) + ' FOO BAR ' + random.choice(['1', '2']))

    # Individual string
    ind_string = 'present.' + str(i + 1)

    # Build both filepaths
    chr_path_1 = pre_dir + ind_string + '.1.fa'
    chr_path_2 = pre_dir + ind_string + '.2.fa'

    # Open both chromosome files
    with open(chr_path_1, 'r') as chr1, open(chr_path_2, 'r') as chr2:
        # Extract sequences
        haplotypes.append(chr1.readlines()[1])
        haplotypes.append(chr2.readlines()[1])

# Get previously computed list with all segsites
all_segsites = []
with open(segsite_filepath, 'r') as segsite_file:
    for line in segsite_file:
        all_segsites.append(int(line))

# Cleanup
os.remove(segsite_filepath)

# List with SNPs (index, [alleles])
snps = []
for site in all_segsites:
    ref_site = haplotypes[0][site]

    if ref_site == '\n':
        continue

    alleles = [ref_site]
    snp = (site, alleles)

    for hap in haplotypes:
        if not hap[site] in alleles:
            alleles.append(hap[site])

    # Ignore triallelic+ sites
    if len(snp[1]) == 1:
        snp[1].append(ref_site)

    if len(snp[1]) == 2:
        snps.append(snp)

# Begin output creation of legend and hap files

# LEGEND and HAP file
for (snp_index, snp) in zip(range(1, len(snps) + 1), snps):

    # LEGEND file
    leg_line = 'SNP' + str(snp_index) + ' '
    leg_line += str(snp[0] + 1) + ' '
    leg_line += snp[1][0] + ' ' + snp[1][1]
    leg_f.append(leg_line)

    # HAP file
    hap_line = ''

    for hap in haplotypes:
        hap_base = hap[snp[0]]

        if hap_base == snp[1][0]:
            hap_line += '0 '
        elif hap_base == snp[1][1]:
            hap_line += '1 '

    # Remove trailing space
    hap_f.append(hap_line[:-1])

# Write to files
with open(f'{out_dir}ref_panel.sam', 'w') as sam, open(f'{out_dir}ref_panel.leg', 'w') as leg, open(f'{out_dir}ref_panel.hap', 'w') as hap:
    for line in sam_f:
        sam.write(line + '\n')
    for line in leg_f:
        leg.write(line + '\n')
    for line in hap_f:
        hap.write(line + '\n')
