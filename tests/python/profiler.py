#!/usr/env python

"""Run a basic simulation of the example origami system."""

import cProfile
import pstats
import sys
sys.path.append('../')
from lattice_dna_origami.lattice_origami_domains import *

# Specificy initial configuration by setting input file and step number
#input_file = JSONInputFile('single_domain.json')
#input_file = JSONInputFile('two_domain.json')
#input_file = JSONInputFile('four_domain_loop.json')
input_file = JSONInputFile('snodin_unbound.json')
#input_file = JSONInputFile('tile_unbound.json')

#system = sys.argv[1]
#moveset = sys.argv[2]
system = 'snodin'
moveset = 'cbmc'
output_file_name = 'profile.hdf5'

filename_base = '{}_{}'.format(system, moveset)
profile_filename = filename_base + '.stats'
stats_filename = '{}.txt'.format(filename_base)

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
origami_system = OrigamiSystemSixteen(input_file, step, temp, strand_M, cation_M)

# Specify moves to be used and associated probabilities
mmc_moveset = {MOVETYPE.EXCHANGE_STAPLE: 0.25,
               MOVETYPE.REGROW_STAPLE: 0.25,
               MOVETYPE.REGROW_SCAFFOLD: 0.25,
               MOVETYPE.ROTATE_ORIENTATION_VECTOR: 0.25}

cbmc_moveset = {MOVETYPE.EXCHANGE_STAPLE: 0.25,
                MOVETYPE.CB_REGROW_STAPLE: 0.25,
                MOVETYPE.CB_REGROW_SCAFFOLD: 0.25,
                MOVETYPE.ROTATE_ORIENTATION_VECTOR: 0.25}


ctcbmc_moveset = {MOVETYPE.EXCHANGE_STAPLE: 0.25,
                  MOVETYPE.CB_REGROW_STAPLE: 0.25,
                  MOVETYPE.CB_CONSERVED_TOPOLOGY: 0.25,
                  MOVETYPE.ROTATE_ORIENTATION_VECTOR: 0.25}

move_settings = cbmc_moveset

# Specify output file type and name

output_file = HDF5OutputFile(output_file_name, origami_system,
        config_write_freq=10000,
        count_write_freq=0,
        energy_write_freq=0)

# Setup up simulation
sim = GCMCSimulation(origami_system, move_settings, output_file,
        center_freq=10000)

# Profile
N = 100000
command_string = 'sim.run(N, logging=100)'
cProfile.run(command_string, profile_filename)

# Analyze results
sys.stdout = open(stats_filename, 'w')
profile_stats = pstats.Stats(profile_filename)
profile_stats.strip_dirs()
profile_stats.sort_stats('time')
profile_stats.print_stats()
