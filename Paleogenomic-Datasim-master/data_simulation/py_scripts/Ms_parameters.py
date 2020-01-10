#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 07:09:10 2019

@author: Liliana
"""

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



import msprime, sys

FileToOpen = sys.argv[1]

with open(FileToOpen) as f:
    data ={}
    for line in f:
        key,value = line.strip().split("=")
        data[key] = float(value)
print(data)

Ne = int(data["Ne"])
length = int(data["length"])
recombination = data["recombination"]
mutation_rate = data["mutation_rate"]
replicates = int(data["replicates"])
modern_samples = int(data["modern_samples"])
ancient_samples = int(data["ancient_samples"])
time = int(data["time"])
num_bases = 10000
tree_filepath = 'tree_data'
newick_filepath = 'newick_tree'


print(Ne)
print(length)
print(recombination)
print(mutation_rate)
print(replicates)
print(modern_samples)
print(ancient_samples)
print(time)

samples = []
if modern_samples > 0:
        for x in range(0,modern_samples):
            samples.append(msprime.Sample(0,0))
if ancient_samples > 0:
        for y in range(0,ancient_samples):
            samples.append(msprime.Sample(0,time))
print(samples)

tree_seq = msprime.simulate(Ne=Ne, mutation_rate=mutation_rate, samples=samples, length=length)

tree = tree_seq.first()
print(tree.draw(format="unicode"))

afs = tree_seq.allele_frequency_spectrum(polarised=True, span_normalise=False)
print(afs)

partitions = write_newick(tree_seq, num_bases, newick_filepath, tree_filepath)
