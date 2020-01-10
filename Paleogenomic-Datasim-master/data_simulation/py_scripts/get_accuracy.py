#####################################################
# This script will analyze the accuracy of SHAPEIT  #
# phasing on simulated data                         #
#####################################################

import os

# Working dir
data_dir = os.path.abspath(os.path.join(os.getcwd(), '../../'))
phase_dir = os.path.abspath(os.path.join(os.getcwd(), 'phased/'))

# Make sure directory exists
if not os.path.exists(phase_dir):
    print('phasing data folder not present')
    quit()

# Get file list
files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if
        os.path.isfile(os.path.join(data_dir, f))]

# Make sure it has the correct format, we need two chromosomes
# per individual
print(data_dir)
print(phase_dir)
if data_dir + '/chr.1.fa' not in files or data_dir + '/chr.2.fa' not in files:
    print('Simulated chromosome files not present')
    quit()

# Get chromosome filenames
chr_path_1 = data_dir + '/chr.2.fa'
chr_path_2 = data_dir + '/chr.1.fa'

# Get phasing data
files = [os.path.join(phase_dir, f) for f in os.listdir(phase_dir) if
        os.path.isfile(os.path.join(phase_dir, f))]

if phase_dir + '/readset.phased.haps' not in files:
    print('Phasing results not present')
    quit()

haps_path = phase_dir + '/readset.phased.haps'

# Lines that will be written to statistics file
out = []
out_alt = []
hits = 0
alt_hits = 0

# Open both chromosome files
with open(chr_path_1, 'r') as chr1_f, open(chr_path_2, 'r') as chr2_f, open(haps_path, 'r') as haps_f:
    # Extract sequences
    chr1 = chr1_f.readlines()[1]
    chr2 = chr2_f.readlines()[1]

    # Extract phasing results
    haps = [x.split() for x in haps_f.readlines()]
    haps = [[x[2], x[3], x[4], x[5], x[6]] for x in haps]

for hap in haps:
    site = int(hap[0])
    allele1 = hap[1]
    allele2 = hap[2]
    hap1_res = int(hap[3])
    hap2_res = int(hap[4])

    # Get site at actual haplotypes
    hap1_site = chr1[site - 1]
    hap2_site = chr2[site - 1]

    # Calculate expected results
    expected1 = -1
    if hap1_site == allele1:
        expected1 = 0
    elif hap1_site == allele2:
        expected1 = 1

    expected2 = -1
    if hap2_site == allele1:
        expected2 = 0
    elif hap2_site == allele2:
        expected2 = 1

    if expected1 == -1 or expected2 == -1:
        print('Triallelic site detected at pos:' + str(site))
        continue

    # Compare to SHAPEIT results
    hit1 = expected1 == hap1_res
    hit2 = expected2 == hap2_res
    hit1_alt = expected1 == hap2_res
    hit2_alt = expected2 == hap1_res

    if hit1:
        hits += 1
    if hit2:
        hits += 1
    if hit1_alt:
        alt_hits += 1
    if hit2_alt:
        alt_hits += 1

    out.append('SHAPEIT: %d, %s, %s, %d, %d | EXP: %d, %d | %s, %s\n' % (site, allele1, allele2, hap1_res, hap2_res, expected1, expected2, hit1, hit2))
    out_alt.append('SHAPEIT: %d, %s, %s, %d, %d | EXP: %d, %d | %s, %s\n' % (site, allele1, allele2, hap1_res, hap2_res, expected2, expected1, hit1_alt, hit2_alt))

out_success = 0
out_alt_success = 0
if len(out) != 0:
    out_success = float(hits) * 100.0 / float(len(out) * 2)
if len(out_alt) != 0:
    out_alt_success = float(alt_hits) * 100.0 / float(len(out_alt) * 2)

with open(os.getcwd() + '/phasing_results.txt', 'w') as out_f:
    if out_success > out_alt_success:
        out_f.writelines(out)
        out_f.write('SUCCESS RATE: {} PERCENT\n'.format(out_success))
    else:
        out_f.writelines(out_alt)
        out_f.write('SUCCESS RATE: {} PERCENT\n'.format(out_alt_success))
