#!/usr/bin/env python

import argparse
import pdb
import sys
#sys.path.insert(0, '../../src/lattice_origami_domains')
sys.path.insert(0, '../lattice_origami_domains')

from lattice_dna_origami.origami_io import *

class OrigamiShell:
    pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('treps', type=int, help='Number of temperature reps')
    parser.add_argument('writes', type=int, help='Number of configurations in traj files')
    parser.add_argument('system_filename', type=str, help='System file')
    parser.add_argument('filebase', type=str, help='Base name for files')

    args = parser.parse_args()
    writes = args.writes
    treps = args.treps
    system_filename = args.system_filename
    filebase = args.filebase

    # Setup files
    # Keys are rep number
    traj_files = {}

    # Swap file
    swap_filename = filebase + '.swp'
    with open(swap_filename) as inp:
        swapsraw = inp.read().splitlines()

    temps = [int(temp) for temp in swapsraw[0].split()]
    swap_data = np.loadtxt(swap_filename, skiprows=1)
    swaps = {}
    for col, temp in enumerate(temps):
        swaps[temp] = swap_data[:, col]

    # Keys are temperatures
    output_files = {}
    for rep in range(treps):
        traj_filename = '{}-{}.trj'.format(filebase, rep)
        input_file = JSONInputFile(system_filename)
        traj_file = PlainTextTrajFile(traj_filename, input_file)
        traj_files[rep] = traj_file

        temp = temps[rep]
        output_filename = '{}-{}.trj'.format(filebase, int(temp))
        output_file = PlainTextTrajOutFile(output_filename)
        output_files[temp] = output_file

    num_swaps = len(swaps[temp])
    swap_freq = writes // num_swaps
    swap = 0
    for write in range(writes):
        for temp in temps:
            rep = swaps[temp][swap]
            traj_file = traj_files[rep]
            output_file = output_files[temp]

            output_file.write_config(traj_file.chains(write), write)

        if (write + 1) % swap_freq == 0:
            swap += 1

if __name__ == '__main__':
    main()
