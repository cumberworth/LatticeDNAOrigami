#!/usr/bin/env python

"""Extract final configurations for restarting REMC runs."""

import argparse
import pdb
import sys

from origamipy import files


def main():
    args = parse_args()
    swap_inp_filename = '{}.swp'.format(args.inp_filebase)
    swap_file = files.UnparsedSingleLineStepInpFile(swap_inp_filename, headerlines=1)
    swap_out_filename = '{}.swp.restart'.format(args.out_filebase)
    with open(swap_out_filename, 'w') as out:
        out.write(swap_file.header)
        out.write(swap_file.get_last_step())

    for thread in range(args.threads):
        traj_inp_filename = '{}-{}.trj'.format(args.inp_filebase, thread)
        traj_file = files.UnparsedMultiLineStepInpFile(traj_inp_filename)
        traj_out_filename = '{}-{}.trj.restart'.format(args.out_filebase, thread)
        with open(traj_out_filename, 'w') as out:
            out.write(traj_file.get_last_step())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'inp_filebase',
        type=str,
        help='Input filebase')
    parser.add_argument(
        'out_filebase',
        type=str,
        help='Output filebase')
    parser.add_argument(
        'threads',
        type=int,
        help='Number of threads/replicas')

    return parser.parse_args()


if __name__ == '__main__':
    main()
