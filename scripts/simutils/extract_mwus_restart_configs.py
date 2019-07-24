#!/usr/bin/env python3

"""Extract specified file of given trajectory to single config trajfile"""

import argparse
import pdb
import sys

from origamipy import files
from origamipy import us_process

def main():
    args = parse_args()
    tags, wins = us_process.read_windows_file(args.wins_filename)
    for i, win in enumerate(wins):
        inp_postfix = '_iter-{}.trj'.format(args.iteration)
        trj_inp_filename = us_process.create_win_filename(win,
                args.inp_filebase, inp_postfix)
        trj_out_filename = us_process.create_win_filename(win,
                args.out_filebase, '.trj.restart')
        trj_file = files.UnparsedMultiLineStepInpFile(trj_inp_filename)
        with open(trj_out_filename, 'w') as out:
            out.write(trj_file.get_last_step())


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
            'wins_filename',
            type=str,
            help='Windows file')
    parser.add_argument(
            'iteration',
            type=str,
            help='Iteration')

    return parser.parse_args()


if __name__ == '__main__':
    main()
