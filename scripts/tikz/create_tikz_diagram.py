#!/usr/bin/env python

"""Generates tikz scripts for configurations of the origami model.

Takes template tikz scripts for geneting configuration diagrams and outputs a
script for configurations at specified steps in the specified intput file.
"""

import argparse
import sys
import string

#sys.path.insert(0, '../../../../lib/lattice_origami_domains')
sys.path.insert(0, '../../../lattice_origami_domains')

from lattice_dna_origami.origami_io import JSONInputFile, PlainTextTrajFile


TEMPLATE_FILENAME = '../../../../share/tikz_template.tex'


def make_tikz_position_bond_orientation_list(chains):
    scaffold_list = ''
    staple_list = ''
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

            bix = chain['orientations'][domain_index][0] * 0.5
            biy = chain['orientations'][domain_index][1] * 0.5
            biz = chain['orientations'][domain_index][2] * 0.5

            if chain['identity'] == 0:
                scaffold_list = scaffold_list + '{} / {} / {} / {} / {} / {} / {} / {} / {}, '.format(
                        rix, riy, riz, aix, aiy, aiz, bix, biy, biz)
            else:
                staple_list = staple_list + '{} / {} / {} / {} / {} / {} / {} / {} / {}, '.format(
                        rix, riy, riz, aix, aiy, aiz, bix, biy, biz)

    # Remove last comma and space
    scaffold_list = scaffold_list[:-2]
    staple_list = staple_list[:-2]

    return scaffold_list, staple_list


def insert_list_and_write(scaffold_list, staple_list, output_filename, coor1,
        coor2, axis):
    with open(TEMPLATE_FILENAME) as input:
        template = string.Template(input.read())

    template = template.substitute(scaffold_list=scaffold_list,
            staple_list=staple_list, coor1=coor1, coor2=coor2, axis=axis)
    with open(output_filename, 'w') as output:
        output.write(template)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('system_filename',
            help='System filename')
    parser.add_argument('traj_filename',
            help='Trajectory filename')
    parser.add_argument('output_filename',
            help='Tex filename')
    parser.add_argument('step', type=int,
            help='Step in configuration file to draw')
    parser.add_argument('--coor1', dest='coor1', default=40,
            help='First perspective coordinate (default 40)')
    parser.add_argument('--coor2', dest='coor2', default=110,
            help='First perspective coordinate (default 110)')
    parser.add_argument('--noaxis', dest='noaxis', action='store_true',
            default=False,
            help='Switch off axis')

    args = parser.parse_args()
    system_filename = args.system_filename
    traj_filename = args.traj_filename
    output_filename = args.output_filename
    step = args.step
    coor1 = args.coor1
    coor2 = args.coor2
    noaxis = args.noaxis
    if noaxis:
        axis = '%'
    else:
        axis = ''
        

    system_file = JSONInputFile(system_filename)
    traj_file = PlainTextTrajFile(traj_filename, system_file)
    chains = traj_file.chains(step)
    scaffold_list, staple_list = make_tikz_position_bond_orientation_list(chains)
    insert_list_and_write(scaffold_list, staple_list, output_filename, coor1,
            coor2, axis)

if __name__ == '__main__':
    main()
