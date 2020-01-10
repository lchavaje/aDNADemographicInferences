#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 07:09:10 2019

@author: Liliana
"""

with open("Msprime_parameters_1.txt") as f:
    data ={}
    for line in f:
        key,value = line.split("=")
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

print(Ne)
print(length)
print(recombination)
print(mutation_rate)
print(replicates)
print(modern_samples)
print(ancient_samples)
print(time)

