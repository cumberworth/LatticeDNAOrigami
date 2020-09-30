#!/usr/bin/env python3

"""Select starting configs for each US window"""

import argparse
import copy
import json
import random
import sys

from origamipy import datatypes
from origamipy import files
from origamipy import us_process


def main():
    args = parse_args()

#    random.seed(123450897)
    traj_file = files.TxtTrajInpFile(
        args.sim_filebase + '.trj', args.system_file)
    bias_tags, wins = us_process.read_windows_file(args.win_file)
    bias_functions = json.load(open(args.bias_functions_filename))
    op_tags = us_process.get_op_tags_from_bias_functions(
        bias_functions, bias_tags)
    ops = datatypes.OrderParams.from_file(args.sim_filebase)
    op_to_config = sort_by_ops(ops, op_tags)
    for i, win in enumerate(wins):
        config = us_process.select_config_by_op_in_win(
            i, win, ops, op_to_config)
        win_filename = us_process.create_win_filename(
            win, args.out_filebase, args.ext)
        config_file = files.TxtTrajOutFile(win_filename)
        config_file.write_config(traj_file.get_chains(config), 0)


def sort_by_ops(ops, tags):
    """Return dictionary of ops to original indices"""
    op_to_config = {}
    for i in range(ops.steps):
        op_key = tuple([ops[tag][i] for tag in tags])
        if op_key in op_to_config:
            op_to_config[op_key].append(i)
        else:
            op_to_config[op_key] = [i]

    return op_to_config


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'system_file',
        type=str,
        help='System file')
    parser.add_argument(
        'win_file',
        type=str,
        help='Windows file')
    parser.add_argument(
        'bias_functions_filename',
        type=str,
        help='Bias functions filename')
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

    return parser.parse_args()


if __name__ == '__main__':
    main()
