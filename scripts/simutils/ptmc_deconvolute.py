#!/usr/bin/env python

import argparse
import pdb
import sys
sys.path.insert(0, '../../src/lattice_origami_domains')

from lattice_dna_origami.origami_io import *

class OrigamiShell:
    pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('treps', type=int, help='Number of temperature reps')
    parser.add_argument('writes', type=int, help='Number of configurations in traj files')
    parser.add_argument('swap_burn_in', type=int, help='Number of steps to discard after swap')
    parser.add_argument('system_filename', type=str, help='System file')
    parser.add_argument('filebase', type=str, help='Base name for files')

    args = parser.parse_args()
    writes = args.writes
    swap_burn_in = args.swap_burn_in
    treps = args.treps
    system_filename = args.system_filename
    filebase = args.filebase

    # Setup files
    # Keys are rep number
    traj_files = {}
    counts_files = {}
    out_files = {}
    op_files = {}

    # Swap file
    swap_filename = filebase + '.swp'
    with open(swap_filename) as inp:
        swapsraw = inp.read().splitlines()

    #temps = [int(temp) for temp in swapsraw[0].split()]
    temps = [330,335,340,345,350,355,360]
    swap_data = np.loadtxt(swap_filename, skiprows=1)
    swaps = {}
    for col, temp in enumerate(temps):
        swaps[temp] = swap_data[:, col]

    # Keys are temperatures
    output_files = {}
    counts_output_files = {}
    out_output_files = {}
    op_output_files = {}
    for rep in range(treps):
        traj_filename = '{}-{}.trj'.format(filebase, rep)
        counts_filename = '{}-{}.counts'.format(filebase, rep)
        out_filename = '{}-{}.out'.format(filebase, rep)
        ops_filename = '{}-{}.order_params'.format(filebase, rep)
        counts = []
        out = []
        ops = []
        with open(counts_filename) as inp:
            counts.append(inp.readlines())

        with open(out_filename) as inp:
            out.append(inp.readlines())

        with open(ops_filename) as inp:
            ops.append(inp.readlines()[1:])

        input_file = JSONInputFile(system_filename)
        traj_file = PlainTextTrajFile(traj_filename, input_file)
        traj_files[rep] = traj_file
        counts_files[rep] = counts
        out_files[rep] = out
        op_files[rep] = ops

        temp = temps[rep]
        output_filename = '{}-{}.trj'.format(filebase, int(temp))
        counts_output_filename = '{}-{}.counts'.format(filebase, int(temp))
        out_output_filename = '{}-{}.out'.format(filebase, int(temp))
        op_output_filename = '{}-{}.order_params'.format(filebase, int(temp))
        output_file = PlainTextTrajOutFile(output_filename)
        output_files[temp] = output_file
        counts_output_files[temp] = open(counts_output_filename, 'w')
        out_output_files[temp] = open(out_output_filename, 'w')
        op_output_files[temp] = open(op_output_filename, 'w')

    num_swaps = len(swaps[temp])
    if num_swaps - 1 > writes:
        swap_inc = (num_swaps - 1) // writes
        swap_freq = 1
    else:
        swap_inc = 1
        swap_freq = writes // (num_swaps - 1)
    swap = swap_inc - 1
    write = swap_burn_in
    while write < writes:
        for temp in temps:
            rep = swaps[temp][swap]
            traj_file = traj_files[rep]
            counts_file = counts_files[rep]
            out_file = out_files[rep]
            op_file = op_files[rep]
            output_file = output_files[temp]
            counts_output_file = counts_output_files[temp]
            out_output_file = out_output_files[temp]
            op_output_file = op_output_files[temp]

            output_file.write_config(traj_file.chains(write), write)
            counts_output_file.write(counts_file[0][write])
            out_output_file.write(out_file[0][write])
            for i in range(10):
                op_output_file.write(op_file[0][write*10 + i])

        if (write + 1) % swap_freq == 0:
            swap += swap_inc
            write += swap_burn_in

        write += 1

if __name__ == '__main__':
    main()
