import msprime
import sys
import subprocess
import shutil
import os
import math
#import argparse


import parameters_input as pt #Must be in the same directory as this file
print('\n*********************')
print(' Input parameters')
print('*********************\n')
sample_size = pt.sample_size
print("sample size = ",sample_size)
Ne = pt.Ne
print("Ne = ",pt.Ne)
length = pt.length
print("length = ",length)
recombination_rate = pt.recombination_rate
print("recombination_rate = ",recombination_rate)
recombination_map = pt.recombination_map
mutation_rate = pt.mutation_rate
print("mutation_rate = ",mutation_rate)
#migration_matrix = pt.migration_matrix
#print(migration_matrix)
demographic_events = pt.demographic_events
print("demographic_events = ",demographic_events)
num_replicates = pt.num_replicates
print("num_replicates = ",num_replicates)
modern_samples = pt.modern_samples
print("This is the number of modern samples plus contaminant and reference = ",modern_samples)
ancient_samples = pt.ancient_samples
print("This is the number of ancient samples = ",ancient_samples)
time = pt.time
print("time =",time)
branch_scale=pt.branch_scale
print("branch_scale = ",branch_scale)
total_samples = 2*(ancient_samples + modern_samples)
print("total samples =",total_samples)

# Generate coalescent trees with msprime
def msprime_sim(sample_size, Ne, length, recombination_rate, 
                recombination_map, mutation_rate,
                demographic_events, num_replicates):

     # Create sample list
    samples = []
    for x in range(modern_samples):
            samples.append(msprime.Sample(0,0))
            samples.append(msprime.Sample(0,0))
    for z in range(0,ancient_samples):
            samples.append(msprime.Sample(0,time))
            samples.append(msprime.Sample(0,time))   
    #print("samples = ",samples)

    # Create population_configurations list
    population_configurations = []
    population_configurations.append(msprime.PopulationConfiguration())
    print("population_configurations =",population_configurations)
    print('\n*******************************')
    print('\n')
    # Run simulation and extract results
    tree_seq = msprime.simulate(sample_size=sample_size, Ne=Ne, 
                            length=length,
                            recombination_rate=recombination_rate,
                            recombination_map=recombination_map, 
                            mutation_rate=mutation_rate, 
                            population_configurations=
                            population_configurations, 
                            demographic_events=demographic_events,
                            samples=samples, 
                            num_replicates=num_replicates)
    
    return tree_seq

#Print jfs
def get_jfs(tree_seq):

    present_nodes = []
    ancient_nodes = []
    x = 2
    y = 64

    while int(x) < 62 :
        present_nodes.append(x)
        x = int(x) + 1

    while int(y) < 74  :
        ancient_nodes.append(y)
        y = int(y) + 1
    print(present_nodes)
    print(ancient_nodes)
    jfs = tree_seq.allele_frequency_spectrum([present_nodes,ancient_nodes],
                                             polarised=True,
                                             span_normalise=False)

    return(jfs)


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
    #print("These are the resulting intervals:", intervals)
    print('\n')
    # Fix rounding error
    diff = length - sum(intervals)

    if diff != 0:
        intervals[len(intervals) - 1] += diff

    #print("These are the intervals with the error fixed:", intervals)

    # Get number of partitions and add intervals
    partitions = 0
    added_intervals = []
    with open(newick_filepath, 'r') as newick_file:
        for line, interval in zip(newick_file, intervals):
            added_intervals.append('[' + str(interval) + '] ' + line)
            partitions += 1
    #print("Added intervals:", added_intervals)

    # Overwrite Newick file with added intervals
    with open(newick_filepath, 'w') as newick_file:
        newick_file.writelines(added_intervals)

    print("Number of partitions =", partitions)
    print("Intervals =", intervals)
    return partitions, intervals

#Split newick_file according to partitions

def split_newick(newick_filepath,newick_dir):
    # Divide newick_file by partitions
    partition_number = 1
    with open(newick_filepath, 'r') as newick_file:
        for line in newick_file:
         line = line.strip(";")
	 # Write files for each partition
         filename ='newicktree_'+str(partition_number)
         partition_filepath = f'{newick_dir}/'+filename
         with open(partition_filepath, 'w') as f:
             f.write(line)
         partition_number +=1    
            
            
# Run seq-gen on Newick tree
def seqgen_sim(partitions,intervals, newick_dir, seq_dir, branch_scale, total_samples, all_sequences_filepath):            

    for i in range(partitions):

        # Define filepaths and variables
        sequence_file = 'sequence_data_'+str(i+1)
        seqgen_filepath = seq_dir+sequence_file
        newicktree_file = 'newicktree_'+str(i+1)
        newick_tree_filepath =newick_dir+newicktree_file
        aux_1 = 1
        p_length = intervals[i]
   
        #run seq-gen
        with open(seqgen_filepath, 'w') as seqgen_file:
            subprocess.run(['seq-gen', '-mHKY', '-l' + str(p_length),
            '-s' + str(branch_scale), '-wa', '-p', str(aux_1),
            newick_tree_filepath], stdout=seqgen_file)
    

    # Sort sequences, msprime does not output chromosomes in order.
    # We will also remove the header, since this sequence file will be
    #split into several files representing our individuals
    
    chr_sequences = []
    master_list = []
    aux_list = []    
    
    for k in range(partitions):
        #Read each sequence file and insert in list
        seqdata_file = 'sequence_data_'+str(k+1)
        seqdata_filepath = seq_dir+seqdata_file 
        with open(seqdata_filepath, 'r') as file:
            chr_sequences = file.readlines()

    	 # Remove header and sort
        chr_sequences.pop(0)
        chr_sequences.sort(key=lambda s : int(s.split()[0]))
         
        # Remove unwanted sequences
        del chr_sequences[total_samples:len(chr_sequences)-1]
         
	 # Write only sequences to file
        with open(seqdata_filepath, 'w') as seqgen_file:
            for line in chr_sequences:
                seqgen_file.write(line.split()[1] + '\n')

        with open(seqdata_filepath, 'r') as seqgen_file:
            chr_sequences = seqgen_file.readlines() 
       
	 #Remove enters from all the elements except for the last partition
        if k < (partitions - 1):
            j = 0
            for line in chr_sequences:
                chr_sequences[j]=line.strip()
                j += 1
	
	 #Create a master list with all the partitions
        if k == 0:
            master_list = chr_sequences[:]
            aux_list = chr_sequences[:]
           # print(master_list)
           # print(aux_list)
        else:
            for i in range (len(chr_sequences)):
                aux_list[i]= master_list[i] + chr_sequences[i]
                del master_list[:]
                master_list = aux_list[:]
           
       # print(master_list)
        del chr_sequences[:]
       
        #Write all the sequences to a file
        with open(all_sequences_filepath, 'w') as file:
            for element in master_list:
                file.write(element)

def split_sequences(all_sequences_filepath, segsite_filepath, modern_samples, ancient_samples, length, dirs):
    # Extract dirs
    con_dir = dirs[0]
    pre_dir = dirs[1]
    anc_dir = dirs[2]
    ref_dir = dirs[3]
    sta_dir = dirs[4]
 

    # Get segsites in simulated sequences, for diagnostics
    # and reference panel
    seqgen_segsites = 0
    with open(all_sequences_filepath, 'r') as seqgen_file, open(segsite_filepath, 'w') as segsite_file:
        chr_sequences = seqgen_file.readlines()

        for i in range(length):
            start = chr_sequences[0][i]
            for line in chr_sequences[1:]:

                # Found a segsite
                if line[i] != start:
                    seqgen_segsites += 1
                    segsite_file.write(f'{i + 1}\n')
                    break

    #print(chr_sequences)

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
    for i in range(modern_samples-1):
        # Individual string
        ind_string = 'present.' + str(i + 1)
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


    #############################
    # Split extra sequences     # 
    #############################
 
    # Individual string
    filename = 'ancestral.fa'
    with open(sta_dir + filename, 'w') as f:
    	f.write('>ref_1' + '\n')
    	f.write(chr_sequences[chr_index])

    return seqgen_segsites

def main():
    
    
    ########################
    # Get coalescent trees #
    ######################## 
       tree_seq = msprime_sim(sample_size, Ne, length, recombination_rate, 
                recombination_map, mutation_rate, demographic_events, num_replicates)
        
       msprime_segsites = tree_seq.get_num_mutations()
   
       # print(tree_seq.tables.nodes)  
  
       #Print joint frequency spectrum
       jfs = get_jfs(tree_seq)
       print("Joint frequency spectrum :")
       print('\n')
       print(jfs)
       print('\n')
       print('\n***************************')

       g_dir = sys.argv[2]
       # print ( "My variable g_dir is equal to ", g_dir )
    
      ##############################
      # Run seq-gen on Newick tree #
      ##############################
       tree_filepath = f'{g_dir}/tree_data'
       newick_filepath = f'{g_dir}/newick_tree'
       segsite_filepath = f'{g_dir}/all_segsites'
       all_sequences_filepath = f'{g_dir}/all_sequences'

       # Transform to Newick tree format
       partitions, intervals = write_newick(tree_seq, length, newick_filepath, 
                                tree_filepath)

       #Create directories for trees & sequences
       newick_dir = f'{g_dir}/newick_trees/'
       seq_dir = f'{g_dir}/sequences/'

       #Create clean directories
       if os.path.exists(newick_dir):
           shutil.rmtree(newick_dir)
       if os.path.exists(seq_dir):
           shutil.rmtree(seq_dir)

       os.makedirs(newick_dir)
       os.makedirs(seq_dir)

       # Split partitions into different files
       split_newick(newick_filepath,newick_dir)

       # Get FASTAs from seq-gen
       seqgen_sim(partitions, intervals, newick_dir, seq_dir, branch_scale, total_samples, all_sequences_filepath) 
          

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
       sta_dir = f'{g_dir}/ancestral_state/'
       dirs = [con_dir, pre_dir, anc_dir, ref_dir, sta_dir]

       # Create clean directories
       for directory in dirs:
           if os.path.exists(directory):
            shutil.rmtree(directory)
           os.makedirs(directory)

        # Put sequences in their respective directories
       seqgen_segsites = split_sequences(all_sequences_filepath, segsite_filepath,modern_samples,
               ancient_samples, length, dirs)
    
        # Cleanup
       if os.path.exists(all_sequences_filepath):
           os.remove(all_sequences_filepath)
       if os.path.exists(newick_dir):
           shutil.rmtree(newick_dir)
       if os.path.exists(seq_dir):
           shutil.rmtree(seq_dir)               

       #####################
       # Output statistics #
       #####################

       # Get a branch scale factor that will result in seq-gen matching msprime
       if branch_scale and seqgen_segsites:
           rec_branch = msprime_segsites / (seqgen_segsites / float(branch_scale))
       else:
           rec_branch = 0

       print('\n')
       print('\n***************************')
       print('input_seqgen.py statistics:')
       print('***************************\n')
       print(f'seq-gen segregating sites = {seqgen_segsites}')
       print(f'msprime segregating sites = {msprime_segsites}')
       print(f'branch scale factor used = {branch_scale}')
       print(f'recommended branch scale factor = {rec_branch}')
       print('\n')

if __name__ == '__main__':
    main()

   
