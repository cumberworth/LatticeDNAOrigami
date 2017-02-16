#!/usr/env python

"""Test input file stuff."""

import sys
sys.path.append('../')
from lattice_dna_origami.lattice_origami_domains import *

input_file = HDF5InputFile('single_domain.hdf5')
step = 0

energies = input_file.energy
print(energies)
