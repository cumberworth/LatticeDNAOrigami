#!/usr/bin/env python3

"""Select starting configs for each US window"""

import argparse
import pdb
import sys
import random
import copy

from origamipy.op_process import read_ops_from_file
from origamipy.op_process import sort_by_ops
from origamipy.origami_io import *
from origamipy.us_process import create_win_filename
from origamipy.us_process import read_windows_file
from origamipy.us_process import select_config_by_op_in_win

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'system_file',
            type=str,
            help='System file')
    parser.add_argument(
            'win_file',
            type=str,
            help='Windows file')
    parser.add_argument(
            'sim_filebase',
            type=str,
            help='Filebase of simulation to select configs from')
    parser.add_argument(
            'out_filebase',
            type=str,
            help='Filebase of starting config trj files')
    parser.add_argument(
            'ext',
            type=str,
            help='Starting config trj file extension')

    args = parser.parse_args()
    system_filename = args.system_file
    wins_filename = args.win_file
    sim_filebase = args.sim_filebase
    out_filebase = args.out_filebase
    ext = args.ext

    traj_file = PlainTextTrajFile(sim_filebase + '.trj', system_filename)
    tags, wins = read_windows_file(wins_filename)
    ops = read_ops_from_file(sim_filebase + '.ops', tags, 0)
    op_to_config = sort_by_ops(ops, tags)
    for i, win in enumerate(wins):
        config = select_config_by_op_in_win(i, win, ops, op_to_config)
        win_filename = create_win_filename(win, out_filebase, ext)
        config_file = PlainTextTrajOutFile(win_filename)
        config_file.write_config(traj_file.chains(config), 0)


if __name__ == '__main__':
    main()
