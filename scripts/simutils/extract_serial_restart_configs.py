#!/usr/bin/env python

"""Extract final configurations for restarting serial runs."""

import argparse
import pdb
import sys

from origamipy import files


def main():
    args = parse_args()
    for temp in args.temps:
        traj_inp_filename = '{}-{}.trj'.format(args.inp_filebase, temp)
        traj_file = files.UnparsedMultiLineStepInpFile(traj_inp_filename)
        traj_out_filename = '{}-{}.trj.restart'.format(args.out_filebase, temp)
        with open(traj_out_filename, 'w') as out:
            out.write(traj_file.get_last_step())


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'inp_filebase',
        type=str,
        help='Input filebase')
    parser.add_argument(
        'out_filebase',
        type=str,
        help='Output filebase')
    parser.add_argument(
        '--temps',
        nargs='+',
        type=int,
        help='Temperatures')

    return parser.parse_args()


if __name__ == '__main__':
    main()
