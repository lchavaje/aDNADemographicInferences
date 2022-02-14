import msprime
import sys
import subprocess
import shutil
import os
import math
#import argparse


import Parameters_test as pt #Must be in the same directory as this file

sample_size = pt.sample_size
print(sample_size)
Ne = pt.Ne
print(pt.Ne)
length = pt.length
print(length)
recombination_rate = pt.recombination_rate
print(recombination_rate)
recombination_map = pt.recombination_map
print(recombination_map)
mutation_rate = pt.mutation_rate
print(mutation_rate)
#print(pt.population_configurations)
migration_matrix = pt.migration_matrix
print(migration_matrix)
demographic_events = pt.demographic_events
print(demographic_events)
num_replicates = pt.num_replicates
print(num_replicates)

modern_samples = pt.modern_samples
print(modern_samples)
ancient_samples = pt.ancient_samples
print(ancient_samples)
time = pt.time
print(time)
num_mod_samples = modern_samples[0] + modern_samples[1]- 1
print(num_mod_samples)

branch_scale=pt.branch_scale
print(branch_scale)

# Generate coalescent trees with msprime
def msprime_sim(sample_size, Ne, length, recombination_rate, 
                recombination_map, mutation_rate, migration_matrix,
                demographic_events, num_replicates):

     # Create sample list
    samples = []

    for w in range(0,len(modern_samples)):
        x = modern_samples[w]
        for y in range (0, x):
            samples.append(msprime.Sample(w,0))
            samples.append(msprime.Sample(w,0))
    for z in range(0,ancient_samples):
            samples.append(msprime.Sample(0,time))
            samples.append(msprime.Sample(0,time))   
    print(samples)

      # Create population_configurations list
    population_configurations = []

    for c in range(0,len(modern_samples)):
            population_configurations.append(
                    msprime.PopulationConfiguration())
    print(population_configurations)
    
    # Run simulation and extract results
    tree_seq = msprime.simulate(sample_size=sample_size, Ne=Ne, 
                            length=length,
                            recombination_rate=recombination_rate,
                            recombination_map=recombination_map, 
                            mutation_rate=mutation_rate, 
                            population_configurations=
                            population_configurations, 
                            migration_matrix=migration_matrix, 
                            demographic_events=demographic_events,
                            samples=samples, 
                            num_replicates=num_replicates)
    
    return tree_seq

# Transform msprime output to Newick tree format
def write_newick(tree_seq, length, newick_filepath, tree_filepath):
    # Output msprime's tree to file
    tree_seq.dump(tree_filepath)

    # Get Newick format tree and partitions
    newick_file = open(newick_filepath, 'w')
    with open(newick_filepath, 'w') as newick_file:
        subprocess.run(['msp', 'newick', '--precision', '14', 
                        tree_filepath], stdout=newick_file)

    # Get each tree's interval, this needs to be appended to the 
        # beginning of each Newick tree described in the file. Intervals 
        # are used by seq-gen to merge the multiple trees that result 
        # from recombination
    intervals = []
    for tree in tree_seq.trees():
        t_length = tree.get_length()
        intervals.append(int(t_length))
    print(intervals)

     # Fix rounding error
    diff = length - sum(intervals)

    if diff != 0:
        intervals[len(intervals) - 1] += diff

    # Get number of partitions and add intervals
    partitions = 0
    added_intervals = []
    with open(newick_filepath, 'r') as newick_file:
        for line, interval in zip(newick_file, intervals):
            added_intervals.append('[' + str(interval) + '] ' + line)
            partitions += 1

    # Overwrite Newick file with added intervals
    with open(newick_filepath, 'w') as newick_file:
        newick_file.writelines(added_intervals)
    
    return partitions

# Run seq-gen on Newick tree
def seqgen_sim(length, branch_scale, partitions, seqgen_filepath, 
               newick_filepath):
    # Run seq-gen
    with open(seqgen_filepath, 'w') as seqgen_file:
        subprocess.run(['seq-gen', '-mHKY', '-l' + str(length),
            '-s' + str(branch_scale), '-p', str(partitions),
            newick_filepath], stdout=seqgen_file)

    # Sort sequences, msprime does not output chromosomes in order.
    # We will also remove the header, since this sequence file will be
    #split into several files representing our individuals
    
    chr_sequences = []
    with open(seqgen_filepath, 'r') as seqgen_file:
        chr_sequences = seqgen_file.readlines()

    # Remove header and sort
    chr_sequences.pop(0)
    chr_sequences.sort(key=lambda s : int(s.split()[0]))

    # Write only sequences to file
    with open(seqgen_filepath, 'w') as seqgen_file:
        for line in chr_sequences:
            seqgen_file.write(line.split()[1] + '\n')

def split_sequences(seqgen_filepath, segsite_filepath, num_mod_samples, ancient_samples, length, dirs):
    # Extract dirs
    con_dir = dirs[0]
    pre_dir = dirs[1]
    anc_dir = dirs[2]
    ref_dir = dirs[3]
    

    # Get segsites in simulated sequences, for diagnostics
    # and reference panel
    seqgen_segsites = 0
    with open(seqgen_filepath, 'r') as seqgen_file, open(segsite_filepath, 'w') as segsite_file:
        chr_sequences = seqgen_file.readlines()

        for i in range(length):
            start = chr_sequences[0][i]
            for line in chr_sequences[1:]:

                # Found a segsite
                if line[i] != start:
                    seqgen_segsites += 1
                    segsite_file.write(f'{i + 1}\n')
                    break

    # Split individuals
    chr_index = 0
   
    #####################
    # Split contaminant #
    #####################
    # First contaminant chromosome
    filename = 'cont.1.fa'
    with open(con_dir + filename, 'w') as f:
        # Write header
        f.write('>cont_1' + '\n')
        # Write only the sequence, no chromosome index
        f.write(chr_sequences[chr_index])

    # Move onto next chromosome
    chr_index += 1

    # Second contaminant chromosome
    filename = 'cont.2.fa'
    with open(con_dir + filename, 'w') as f:
        f.write('>cont_1' + '\n')
        f.write(chr_sequences[chr_index])

    chr_index += 1     
  
    #################################
    # Split present day individuals #
    #################################
    for i in range(num_mod_samples):
        # Individual string
        ind_string = 'present.' + str(i + 1)
        # First chromosome
        filename = ind_string + '.1.fa'
        with open(pre_dir + filename, 'w') as f:
            f.write('>present_' + str(i + 1) + '\n')
            f.write(chr_sequences[chr_index])

        chr_index += 1

        # Second chromosome
        filename = ind_string + '.2.fa'
        with open(pre_dir + filename, 'w') as f:
            f.write('>present_' + str(i + 1) + '\n')
            f.write(chr_sequences[chr_index])

        chr_index += 1


    #############################
    # Split ancient individuals # 
    #############################
    for i in range(ancient_samples):
        # Individual string
        ind_string = 'ancient.' + str(i + 1)
        # First chromosome
        filename = ind_string + '.1.fa'
        with open(anc_dir + filename, 'w') as f:
            f.write('>ancient_' + str(i + 1) + '\n')
            f.write(chr_sequences[chr_index])

        chr_index += 1

        # Second chromosome
        filename = ind_string + '.2.fa'
        with open(anc_dir + filename, 'w') as f:
            f.write('>ancient_' + str(i + 1) + '\n')
            f.write(chr_sequences[chr_index])

        chr_index += 1

    return seqgen_segsites

def main():
    
    
    ########################
    # Get coalescent trees #
    ########################
    tree_seq = msprime_sim(sample_size, Ne, length, recombination_rate, 
                recombination_map, mutation_rate, migration_matrix,
                demographic_events, num_replicates)
     
    msprime_segsites = tree_seq.get_num_mutations()
    print(msprime_segsites)
    
    #Print joint frequency spectrum
    jfs = tree_seq.allele_frequency_spectrum([[0,3],[1,2]],
                                             polarised=True, 
                                             span_normalise=False)
    print(jfs)
    g_dir = sys.argv[2]
    print ( "My variable g_dir is equal to ", g_dir )
    
    ##############################
    # Run seq-gen on Newick tree #
    ##############################
    tree_filepath = f'{g_dir}/tree_data'
    newick_filepath = f'{g_dir}/newick_tree'
    seqgen_filepath = f'{g_dir}/sequence_data'
    segsite_filepath = f'{g_dir}/all_segsites'

    # Transform to Newick tree format
    partitions = write_newick(tree_seq, length, newick_filepath, 
                              tree_filepath)
     # Get FASTAs from seq-gen
    seqgen_sim(length, branch_scale, partitions, seqgen_filepath, 
               newick_filepath)

    # Cleanup
    if os.path.exists(tree_filepath):
        os.remove(tree_filepath)
    if os.path.exists(newick_filepath):
        os.remove(newick_filepath)

    ##############################################
    # Split seq-gen output into individual files #
    ##############################################

    # Directory constants
    con_dir = f'{g_dir}/contaminant/'
    pre_dir = f'{g_dir}/present/'
    anc_dir = f'{g_dir}/ancient/'
    ref_dir = f'{g_dir}/reference/'
    dirs = [con_dir,pre_dir, anc_dir, ref_dir]

    # Create clean directories
    for directory in dirs:
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)

    # Put sequences in their respective directories
    seqgen_segsites = split_sequences(seqgen_filepath, segsite_filepath,num_mod_samples,
            ancient_samples, length, dirs)
    
    # Cleanup
    if os.path.exists(seqgen_filepath):
        os.remove(seqgen_filepath)

    #####################
    # Output statistics #
    #####################

    # Get a branch scale factor that will result in seq-gen matching msprime
    if branch_scale and seqgen_segsites:
        rec_branch = msprime_segsites / (seqgen_segsites / float(branch_scale))
    else:
        rec_branch = 0

    print('\n**************************')
    print('Ms_param_test.py statistics:')
    print('**************************\n')
    print(f'seq-gen segregating sites = {seqgen_segsites}')
    print(f'msprime segregating sites = {msprime_segsites}')
    print(f'branch scale factor used = {branch_scale}')
    print(f'recommended branch scale factor = {rec_branch}')

if __name__ == '__main__':
    main()

   
