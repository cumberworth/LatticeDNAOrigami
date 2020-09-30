#!/usr/bin/env python3

"""Extract specified file of given trajectory to single config trajfile"""

import argparse
import pdb
import sys

from origamipy.origami_io import *
from origamipy.us_process import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inp_filebase', type=str, help='Input filebase')
    parser.add_argument('out_filebase', type=str, help='Output filebase')
    parser.add_argument('wins_file', type=str, help='Windows file')
    parser.add_argument('config', type=int, help='Number of config in file')
    parser.add_argument('iteration', type=str, help='Iteration')
    args = parser.parse_args()
    inp_filebase = args.inp_filebase
    out_filebase = args.out_filebase
    wins_filename = args.wins_file
    step = args.config
    iteration = args.iteration

    tags, wins = read_windows_file(wins_filename)
    for i, win in enumerate(wins):
        inp_postfix = '_iter-{}.trj'.format(iteration)
        inp_win_filename = create_win_filename(win, inp_filebase, inp_postfix)
        out_win_filename = create_win_filename(win, out_filebase, '.restart')
        traj_file = PlainTextTrajFile(inp_win_filename, int(300))
        out_file = PlainTextTrajOutFile(out_win_filename)
        out_file.write_config(traj_file.chains(step), 0)


if __name__ == '__main__':
    main()
