#!/usr/bin/env python

"""Extract config at specified step."""

import argparse
import pdb
import sys

from origamipy import files


def main():
    args = parse_args()
    sys_inp_file = files.JSONStructInpFile(args.sys_filename)
    traj_file = files.TxtTrajInpFile(args.trj_filename, sys_inp_file)
    chains = traj_file.get_chains(args.step)
    sys_out_filename = '{}.json'.format(args.out_filebase)
    sys_out_file = files.JSONStructOutFile(sys_out_filename, sys_inp_file)
    sys_out_file.write(chains)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'sys_filename',
        type=str,
        help='System input filename')
    parser.add_argument(
        'trj_filename',
        type=str,
        help='Trajectory filename')
    parser.add_argument(
        'step',
        type=int,
        help='Step to extract')
    parser.add_argument(
        'out_filebase',
        type=str,
        help='Output filebase')

    return parser.parse_args()


if __name__ == '__main__':
    main()
