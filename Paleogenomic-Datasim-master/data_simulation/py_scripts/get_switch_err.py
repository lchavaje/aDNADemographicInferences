#####################################################
# This script will analyze the accuracy of SHAPEIT  #
# phasing on simulated data                         #
#####################################################

import os

results_path = os.getcwd() + '/phasing_results.txt'

# Make sure directory exists
if not os.path.exists(results_path):
    print('phasing accuracy results not present')
    quit()

# Lines that will be written to statistics file
out = []

# Count switches in results file
switch_counter = 0
phased_sites = 0
switch_flag = False
with open(results_path, 'r') as res_f:
    for line in res_f.readlines():
        if 'SHAPEIT' not in line:
            continue

        phased_sites += 1

        f_block = True if 'False, False' in line else False

        if f_block and switch_flag:
            continue
        elif (f_block and not switch_flag) or (not f_block and switch_flag):
            switch_flag = not switch_flag
            switch_counter += 1

            site = int(line.split()[1].split(',')[0])
            out.append(f'SWITCH AT {site}\n')

if phased_sites == 0:
    switch_rate = 0
else:
    switch_rate = switch_counter / float(phased_sites) * 100.0

with open(os.getcwd() + '/switch_error.txt', 'w') as out_f:
    out_f.writelines(out)
    out_f.write(f'SWITCH RATE: {switch_rate} PERCENT\n')
