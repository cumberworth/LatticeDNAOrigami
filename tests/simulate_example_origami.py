#!/usr/env python

"""Run a basic simulation of the example origami system."""

import sys
sys.path.append('../')
from lattice_dna_origami.lattice_origami_domains import *

# Specificy initial configuration by setting input file and step number
#input_file = JSONInputFile('simple_loop_linear.json')
#input_file = JSONInputFile('simple_loop.json')
#input_file = JSONInputFile('cyclic_example.json')
#input_file = JSONInputFile('single_domain.json')
input_file = JSONInputFile('snodin_unbound.json')
step = 0

# Set conditions
temp = 330

# Staple strand concentration (M)
strand_M = 1e-3

# Cation concentration (M)
#cation_M = 1
cation_M = 0.5

# Setup origami system object
#origami_system = OrigamiSystemEight(input_file, step, temp, cation_M)
origami_system = OrigamiSystemSixteen(input_file, step, temp, strand_M, cation_M)

# Specify moves to be used and associated probabilities
move_settings = {MOVETYPE.CB_EXCHANGE_STAPLE: 0.25,
                 MOVETYPE.CB_REGROW_STAPLE: 0.25,
#                 MOVETYPE.CB_REGROW_SCAFFOLD: 0.25,
                 MOVETYPE.CB_CONSERVED_TOPOLOGY: 0.25,
                 MOVETYPE.ROTATE_ORIENTATION_VECTOR: 0.25}

# Specify output file type and name
#output_file_name = 'simple_loop_replica-0.hdf5'
output_file_name = 'single_domain.hdf5'
config_write_freq = 0
count_write_freq = 1

output_file = HDF5OutputFile(output_file_name, origami_system,
        config_write_freq=config_write_freq,
        count_write_freq=count_write_freq)

# Setup up simulation
sim = GCMCSimulation(origami_system, move_settings, output_file)

# Run
N = 100000
sim.run(N)
