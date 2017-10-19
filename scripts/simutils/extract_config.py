#!/usr/bin/env python3

"""Extract specified config from traj file to new traj file

The main purpose is to create restart traj files that are small for transfer
"""

import argparse
import pdb
import sys

from origamipy.origami_io import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inp_filename', type=str, help='Input filename')
    parser.add_argument('out_filename', type=str, help='Output filename')
    parser.add_argument('step', type=int, help='Step of config in file')

    args = parser.parse_args()

    inp_filename = args.inp_filename
    out_filename = args.out_filename
    step = args.step

    traj_file = PlainTextTrajFile(inp_filename, int(300))
    out_file = PlainTextTrajOutFile(out_filename)
    out_file.write_config(traj_file.chains(step), 0)


if __name__ == '__main__':
    main()
