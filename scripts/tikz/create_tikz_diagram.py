#!/usr/bin/env python

"""Generates tikz scripts for configurations of the origami model.

Takes template tikz scripts for geneting configuration diagrams and outputs a
script for configurations at specified steps in the specified intput file.
"""

import argparse
import string

from origamipy.files import JSONStructInpFile


def main():
    args = parse_args()
    if not args.axis:
        axis = '%'
    else:
        axis = ''

    # This assumes the traj file is json
    traj_file = JSONStructInpFile(args.traj_filename)
    chains = traj_file.chains(args.step)
    inserts = make_tikz_position_bond_orientation_list(
        chains,
        args.orelen,
        args.exclude_scaffold,
        args.exclude_staples)
    insert_list_and_write(
        inserts,
        args.output_filename,
        args.coor1,
        args.coor2,
        axis,
        args.template_filename)


def make_tikz_position_bond_orientation_list(
        chains,
        orelen,
        exclude_scaffold,
        exclude_staples):
    # ndr is for the next domain vector
    scaffold_list = ''
    scaffold_list_ndr = ''
    staple_list = ''
    staple_list_ndr = ''
    for chain in chains:
        for domain_index in range(len(chain['positions'])):
            rix = chain['positions'][domain_index][0]
            riy = chain['positions'][domain_index][1]
            riz = chain['positions'][domain_index][2]
            try:
                rjx = chain['positions'][domain_index + 1][0]
                rjy = chain['positions'][domain_index + 1][1]
                rjz = chain['positions'][domain_index + 1][2]
                aix = rjx - rix
                aiy = rjy - riy
                aiz = rjz - riz
            except IndexError:
                aix = 0
                aiy = 0
                aiz = 0

            bix = chain['orientations'][domain_index][0] * orelen
            biy = chain['orientations'][domain_index][1] * orelen
            biz = chain['orientations'][domain_index][2] * orelen

            entry = '{} / {} / {} / {} / {} / {} / {} / {} / {}, '.format(
                    rix, riy, riz, aix, aiy, aiz, bix, biy, biz)
            if chain['identity'] == 0:
                scaffold_list = scaffold_list + entry
            else:
                staple_list = staple_list + entry

            if domain_index != len(chain['positions']) - 1:
                if chain['identity'] == 0:
                    scaffold_list_ndr = scaffold_list_ndr + entry
                else:
                    staple_list_ndr = staple_list_ndr + entry

    # Remove last comma and space
    scaffold_list = scaffold_list[:-2]
    scaffold_list_ndr = scaffold_list_ndr[:-2]
    staple_list = staple_list[:-2]
    staple_list_ndr = staple_list_ndr[:-2]

    if exclude_scaffold:
        scaffold_list = ''
        scaffold_list_ndr = ''

    if exclude_staples:
        staple_list = ''
        staple_list_ndr = ''

    inserts = {
        'scaffold_list': scaffold_list,
        'scaffold_list_ndr': scaffold_list_ndr,
        'staple_list': staple_list,
        'staple_list_ndr': staple_list_ndr}

    return inserts


def insert_list_and_write(
        inserts,
        output_filename,
        coor1,
        coor2,
        axis,
        template_filename):

    with open(template_filename) as input:
        template = string.Template(input.read())

    template = template.substitute(
        scaffold_list=inserts['scaffold_list'],
        scaffold_list_ndr=inserts['scaffold_list_ndr'],
        staple_list=inserts['staple_list'],
        staple_list_ndr=inserts['staple_list_ndr'],
        coor1=coor1,
        coor2=coor2,
        axis=axis)
    with open(output_filename, 'w') as output:
        output.write(template)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'traj_filename',
        help='Trajectory filename (only json for now)')
    parser.add_argument(
        'output_filename',
        help='Tex filename')
    parser.add_argument(
        'step',
        type=int,
        help='Step in configuration file to draw')
    parser.add_argument(
        '--coor1',
        dest='coor1',
        default=0,
        help='First perspective coordinate (default 0)')
    parser.add_argument(
        '--coor2',
        dest='coor2',
        default=0,
        help='Second perspective coordinate (default 0)')
    parser.add_argument(
        '--orelen',
        dest='orelen',
        default=0.5,
        type=float,
        help='Length of orientation vector (default 0.5)')
    parser.add_argument(
        '--axis',
        dest='axis',
        action='store_true',
        default=False,
        help='Switch off axis')
    parser.add_argument(
        '--exclude_scaffold',
        default=False,
        help='Do not draw scaffold')
    parser.add_argument(
        '--exclude_staples',
        default=False,
        help='Do not draw staples')
    parser.add_argument(
        '--tikz_template',
        dest='template_filename',
        default='tikz_template.tex',
        help='Tikz template')
#    parser.add_argument(
#        '--system_filename',
#        help='System filename if using non-json traj file')

    return parser.parse_args()

if __name__ == '__main__':
    main()
