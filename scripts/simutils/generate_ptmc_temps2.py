#!/usr/bin/env python

"""Generate temperatures along given order parameter."""

import argparse
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.optimize import minimize

from origamipy import files


def main():
    args = parse_args()
    aves = pd.read_csv(args.inp_filename, sep=' ')
    if args.rtag:
        aves = aves[aves[args.rtag] == args.rvalue].reset_index()

    old_temps = aves['temp']
    old_ops = aves[args.tag]

    # Prevent instabilities in minimization (need monotonically decreasing)
    old_ops = old_ops[::-1].sort_values()
    boundary_reached = False
    for i, op in old_ops.items():
        if op > args.max_op*0.95:
            if boundary_reached:
                old_ops = old_ops.drop(i, 0)
                old_temps = old_temps.drop(i, 0)
            else:
                boundary_reached = True

    old_ops = old_ops[::-1]
    boundary_reached = False
    for i, op in old_ops.items():
        if op < args.max_op*0.05:
            if boundary_reached:
                old_ops = old_ops.drop(i, 0)
                old_temps = old_temps.drop(i, 0)
            else:
                boundary_reached = True

#    spline_params = interpolate.splrep(old_temps, old_ops)
    interpolated_ops_f = interpolate.interp1d(old_temps, old_ops, kind='linear',
            fill_value='extrapolate')
    guess_temps = np.linspace(old_temps.iloc[0], old_temps.iloc[-1],
            num=args.threads - 2)
    desired_ops = np.linspace(args.max_op - args.max_op*0.05, args.max_op*0.05, num=args.threads - 2)
    new_temps = minimize(sum_of_squared_errors, guess_temps,
#             args=(desired_ops, spline_params))
             args=(desired_ops, interpolated_ops_f)).x
    new_temps.sort()
    temp_diff = new_temps[-1] - new_temps[0]
    low_temps = [new_temps[0] - temp_diff/4]
    high_temps = [new_temps[-1] + temp_diff/4]
    new_temps = np.concatenate([low_temps, new_temps, high_temps])
    np.set_printoptions(formatter={'float': '{:0.3f}'.format}, linewidth=200)
    new_temps = np.around(new_temps, decimals=3)
    temps_string = ''
    for temp in new_temps:
        temps_string = temps_string + '{:.3f} '.format(temp)

    print(temps_string)


#def sum_of_squared_errors(temps, desired_ops, spline_params):
#    new_ops = interpolate.splev(temps, spline_params, der=0)
#    return ((new_ops - desired_ops)**2).sum()
def sum_of_squared_errors(temps, desired_ops, interpolated_ops_f):
    squared_error = 0
    for temp, op in zip(temps, desired_ops):
        new_op = interpolated_ops_f(temp)
        squared_error += (new_op - op)**2

    return squared_error


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'inp_filename',
        type=str,
        help='Expectations filename')
    parser.add_argument(
        'tag',
        type=str,
        help='Order parameter tag')
    parser.add_argument(
        'max_op',
        type=float,
        help='Maximum value of order parameter')
    parser.add_argument(
        'threads',
        type=int,
        help='Number of threads/replicas')
    parser.add_argument(
            '--rtag',
            type=str,
            help='Tag to slice on')
    parser.add_argument(
            '--rvalue',
            type=float,
            help='Slice value')

    return parser.parse_args()


if __name__ == '__main__':
    main()
