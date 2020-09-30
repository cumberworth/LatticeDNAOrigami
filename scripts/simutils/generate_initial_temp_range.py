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
    temps = np.linspace(args.min_temp, args.max_temp, args.threads)
    temps_string = ''
    for temp in temps:
        temps_string = temps_string + '{:.3f} '.format(temp)

    print(temps_string)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'min_temp',
        type=float,
        help='Minimum temperature')
    parser.add_argument(
        'max_temp',
        type=float,
        help='Maximum temperature')
    parser.add_argument(
        'threads',
        type=int,
        help='Number of threads/replicas')

    return parser.parse_args()


if __name__ == '__main__':
    main()
