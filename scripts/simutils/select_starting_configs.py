#!/usr/bin/env python3

"""Select starting configs for each US window"""

import argparse
import pdb
import sys
import random
import copy
sys.path.insert(0, '../../src/lattice_origami_domains')

from lattice_dna_origami.origami_io import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('system_file', type=str, help='System file')
    parser.add_argument('win_file', type=str, help='Windows file')
    parser.add_argument('sim_filebase', type=str, help='Filebase of simulation to select configs from')
    parser.add_argument('out_filebase', type=str, help='Filebase of starting config trj files')
    args = parser.parse_args()
    system_filename = args.system_file
    wins_filename = args.win_file
    sim_filebase = args.sim_filebase
    out_filebase = args.out_filebase

    traj_file = PlainTextTrajFile(sim_filebase + '.trj', system_filename)
    wins = read_windows(wins_filename)
    ops = parse_counts_file(sim_filebase + '.counts')
    op_to_config = sort_configs_by_ops(ops)
    for i, win in enumerate(wins):
        config = select_config(i, win, ops, op_to_config)
        win_filename = create_win_filename(win, out_filebase)
        config_file = PlainTextTrajOutFile(win_filename)
        config_file.write_config(traj_file.chains(config), 0)


if __name__ == '__main__':
    main()
