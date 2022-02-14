#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 09:48:15 2020

@author: Liliana
"""
#Input data for sample_simulation.py. Both files must be in the same directory.

sample_size=None #(int)If not specified or None, this defaults to the sum of the subpopulation sample sizes. Either sample_size, population_configurations or samples must be specified
Ne=10000 #(float) Defaults to 1 if not specified
length = 100000  #(float) Cannot be used with recombination_map. Dafaults to 1 if not specified
recombination_rate = 1e-8 #(float) Cannot be used with recombination_map. Defaults to 0 if not specified
recombination_map = None #Cannot be used along with recombination_rate or length. Defaults to a uniform rate as described in the recombination_rate
mutation_rate = 1.5e-8 #(float) if not specified, no mutations are generated
#population_configurations = None #(list or None)if not specified, a single population with a sample_size size is assumed
#migration_matrix = [[0,0.00000000001],
#                    [0.000000001,0]] #(list)N x N matrix (from population_configurations N) with 0 on the diagonal
demographic_events = None #(list) Events in non decreasing order
random_seed = None #(int) None generates random seed automatically. Valid random seeds between 1 and 2^32 - 1
num_replicates= None #(int) If not specified or None, no replication is performed
from_ts = None #Check documentation for details. Default is None.
start_time = 0 #(float) Defaults to zero
end_time = None #(float) Default NOne or not specified
record_full_arg = False # (bool) Defaults to False.
model = None # (str or simulation model instance) Check documentation for valid models
#samples (list) Generated in the program. Input data in variables below. Cannot be used in conjunction with sample_size. 
modern_samples = 32 
ancient_samples=5 #Ancient samples is assigned to population 0.
time=20
branch_scale=7e-08 # Branch scaling factor default value=7e-08


