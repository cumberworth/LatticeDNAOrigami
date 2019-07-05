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
    old_temps = aves['temp']
    old_ops = aves[args.tag]

    # Prevent instabilities in minimization (need monotonically decreasing)
    old_ops = old_ops.sort_values()[::-1]
#    spline_params = interpolate.splrep(old_temps, old_ops)
    interpolated_ops_f = interpolate.interp1d(old_temps, old_ops, kind='linear',
            fill_value='extrapolate')
    guess_temps = np.linspace(old_temps[1], old_temps[len(old_temps) - 2],
            num=args.threads - 4)
    desired_ops = np.linspace(args.max_op - 1, 1, num=args.threads - 4)
    new_temps = minimize(sum_of_squared_errors, guess_temps,
#             args=(desired_ops, spline_params))
             args=(desired_ops, interpolated_ops_f)).x
    new_temps.sort()
    low_temps = [new_temps[0] - 3, new_temps[0] - 1, new_temps[0] - 0.3]
    high_temp = new_temps[len(new_temps) - 1]
    high_temps = [high_temp + 1]
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

    return parser.parse_args()


if __name__ == '__main__':
    main()
