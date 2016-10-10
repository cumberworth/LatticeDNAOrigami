#!/usr/env python

"""Run a basic simulation of the example origami system."""

import sys
sys.path.insert(0, '../')

from lattice_dna_origami.origami_io import *
from lattice_dna_origami.origami_system import *


# Specificy initial configuration by setting input file and step number
#input_file = JSONInputFile('simple_loop_linear.json')
#input_file = JSONInputFile('simple_loop.json')
#input_file = JSONInputFile('cyclic_example.json')
#input_file = JSONInputFile('single_domain.json')
#input_file = JSONInputFile('snodin_unbound.json')
input_file = JSONInputFile('snodin_assembled.json')
#input_file = JSONInputFile('four_domain_loop.json')
#input_file = JSONInputFile('two_domain.json')
#input_file = JSONInputFile('two_domain_single.json')
step = 0

# Set conditions
temp = 350

# Staple strand concentration (M)
strand_M = 1e-3

# Cation concentration (M)
#cation_M = 1
cation_M = 1

# Setup origami system object
#origami_system = OrigamiSystemEight(input_file, step, temp, cation_M)
origami_system = OrigamiSystemSixteen(input_file, step, temp, strand_M, cation_M, misbinding=False)

# Specify moves to be used and associated probabilities
move_settings = {MOVETYPE.EXCHANGE_STAPLE: 1}
                 #MOVETYPE.IDENTITY: 1}
                 #MOVETYPE.CB_REGROW_STAPLE: 1}
                 #MOVETYPE.REGROW_SCAFFOLD: 1/4,
                 #MOVETYPE.CB_CONSERVED_TOPOLOGY: 1/4,
                 #MOVETYPE.ROTATE_ORIENTATION_VECTOR: 1/4}

# Specify output file type and name
#output_file_name = 'simple_loop_replica-0.hdf5'
output_file_name = 'test.hdf5'

output_file = HDF5OutputFile(output_file_name, origami_system,
        config_write_freq=10000,
        count_write_freq=0,
        energy_write_freq=0)

# Setup up simulation
sim = GCMCSimulation(origami_system, move_settings, output_file, center_freq=10000)

# Run
N = 100000
sim.run(N, logging=100)
