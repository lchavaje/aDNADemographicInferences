import msprime
import sys
import subprocess
import shutil
import os
import math
import argparse

# Generate coalescent trees with msprime
def msprime_sim(gens, conts, refs, ancs, event, event_time, num_bases):
    # If an event was given, define the necessary PopulationConfigurations
    pop_configs = None
    dem_events = None
    if event == 'split':
        pop_configs = [
            msprime.PopulationConfiguration(initial_size=5e3),
            msprime.PopulationConfiguration(initial_size=5e3)
        ]

        dem_events = [
            msprime.MassMigration(time=event_time, source=0, destination=1, proportion=1.0)
        ]
    elif event == 'bottleneck':
        pop_configs = [
            msprime.PopulationConfiguration(initial_size=1e3),
        ]

        dem_events = [
            msprime.PopulationParametersChange(time=event_time, initial_size=1e4, population_id=0)
        ]

    # Create sample list
    samples = []
    # Add present samples to sample list
    for _ in range(conts + refs):
            samples.append(msprime.Sample(0, 0))
            samples.append(msprime.Sample(0, 0))

    # If there was a population split, we'll use 2 different populations
    if event == 'split':
            ancient_pop = 1
    else:
            ancient_pop = 0

    # Add ancient samples to sample list
    for _ in range(ancs):
            samples.append(msprime.Sample(ancient_pop, gens))
            samples.append(msprime.Sample(ancient_pop, gens))

    # Run simulation and extract results
    if event:
            tree_seq = msprime.simulate(
                    samples=samples, recombination_rate=2e-8,
                    mutation_rate=2e-8, length=num_bases,
                    population_configurations=pop_configs,
                    demographic_events=dem_events)
    else:
            tree_seq = msprime.simulate(
                    samples=samples, recombination_rate=2e-8,
                    mutation_rate=2e-8, length=num_bases, Ne=1e4)

    return tree_seq

# Transform msprime output to Newick tree format
def write_newick(tree_seq, num_bases, newick_filepath, tree_filepath):
    # Output msprime's tree to file
    tree_seq.dump(tree_filepath)

    # Get Newick format tree and partitions
    newick_file = open(newick_filepath, 'w')
    with open(newick_filepath, 'w') as newick_file:
        subprocess.run(['msp', 'newick', '--precision', '14', tree_filepath], stdout=newick_file)

    # Get each tree's interval, this needs to be appended to the beginning
    # of each Newick tree described in the file. Intervals are used by seq-gen
    # to merge the multiple trees that result from recombination
    intervals = []
    for tree in tree_seq.trees():
        length = tree.get_length()
        intervals.append(int(length))

    # Fix rounding error
    diff = num_bases - sum(intervals)

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
def seqgen_sim(num_bases, branch_scale, partitions, seqgen_filepath, newick_filepath):
    # Run seq-gen
    with open(seqgen_filepath, 'w') as seqgen_file:
        subprocess.run(['seq-gen', '-mHKY', '-l' + str(num_bases),
            '-s' + str(branch_scale), '-p', str(partitions),
            newick_filepath], stdout=seqgen_file)

    # Sort sequences, msprime does not output chromosomes in order.
    # We will also remove the header, since this sequence file will be split
    # into several files representing our individuals
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

def split_sequences(seqgen_filepath, segsite_filepath, refs, ancs, num_bases, dirs):
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

        for i in range(num_bases):
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
    for i in range(refs):
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
    for i in range(ancs):
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
    # Important parameters #
    ########################

    # Global dir
    parser = argparse.ArgumentParser(description='Simulate chromosomes that are compatible with gargammel.')

    parser.add_argument('-d', '--dir', nargs=1, help='Directory where the necessary file structures and results will be placed', required=True)
    parser.add_argument('-g', '--gens', nargs=1, help='Age of ancient individuals to be simulated', required=True)
    parser.add_argument('-a', '--ancs', nargs=1, help='Number of ancient individuals to be simulated', required=True)
    parser.add_argument('-r', '--refs', nargs=1, help='Number of present individuals for reference panel', required=True)
    parser.add_argument('-b', '--bases', nargs=1, help='Number of base pairs to be simulated', required=True)
    parser.add_argument('-s', '--branch-scale', nargs=1, help='Branch scaling factor', default=7e-08)
    parser.add_argument('-e', '--event', nargs='?', const=None, choices=['split', 'bottleneck'], help='Type of demographic event to simulate')
    parser.add_argument('-t', '--event-time', nargs='?', const=None, help='Time in generations of the demographic event')

    # Parse the arguments
    args = parser.parse_args()

    # Handle parser errors
    if args.event_time and not args.event:
        parser.error('--event-time requires --event')
    elif args.event and not args.event_time:
        parser.error('--event requires --event-time')

    g_dir = args.dir[0]
    if not os.path.exists(g_dir):
        print(f'{dir} is not a valid directory!')
        return

    gens = int(args.gens[0])
    ancient_individuals = int(args.ancs[0])
    present_individuals = int(args.refs[0])
    num_bases = int(args.bases[0])
    branch_scale = float(args.branch_scale[0])
    if args.event:
        event = args.event
    else:
        event = None
    if args.event_time:
        event_time = int(args.event_time)
    else:
        event_time = None

    # Only 1 contamination individual is currently supported
    contamination_individuals = 1

    ########################
    # Get coalescent trees #
    ########################
    tree_seq = msprime_sim(gens, contamination_individuals, present_individuals,
            ancient_individuals, event, event_time, num_bases)
    msprime_segsites = tree_seq.get_num_mutations()

    ##############################
    # Run seq-gen on Newick tree #
    ##############################
    tree_filepath = f'{g_dir}/tree_data'
    newick_filepath = f'{g_dir}/newick_tree'
    seqgen_filepath = f'{g_dir}/sequence_data'
    segsite_filepath = f'{g_dir}/all_segsites'

    # Transform to Newick tree format
    partitions = write_newick(tree_seq, num_bases, newick_filepath, tree_filepath)
    # Get FASTAs from seq-gen
    seqgen_sim(num_bases, branch_scale, partitions, seqgen_filepath, newick_filepath)

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
    dirs = [con_dir, pre_dir, anc_dir, ref_dir]

    # Create clean directories
    for directory in dirs:
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)

    # Put sequences in their respective directories
    seqgen_segsites = split_sequences(seqgen_filepath, segsite_filepath, present_individuals,
            ancient_individuals, num_bases, dirs)
    
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
    print('gen_genomes.py statistics:')
    print('**************************\n')
    print(f'seq-gen segregating sites = {seqgen_segsites}')
    print(f'msprime segregating sites = {msprime_segsites}')
    print(f'branch scale factor used = {branch_scale}')
    print(f'recommended branch scale factor = {rec_branch}')

if __name__ == '__main__':
    main()
