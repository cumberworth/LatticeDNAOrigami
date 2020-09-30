#!/usr/bin/env python

import json
import pdb
import numpy as np

"""A script to create input json files from sequences."""

STAPLE_SEQFILE = 'tile_staples.seq'
SCAFFOLD_SEQFILE = 'tile_scaffold.seq'
OUTPUT_FILENAME = 'tile_unbound.json'
CYCLIC = True

COMPLEMENTARY_BASE_PAIRS = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def read_seqfile(filename):
    """Read seq files in ??? format and strip metadata."""
    with open(filename) as inp:
        seqs = inp.read().split('\n')
        seqs = [seq for seq in seqs if not '>' in seq]
        seqs = [seq for seq in seqs if seq != '']

    return seqs


def reverse_complement(seq):
    """Return reverse-complement."""
    comp_seq_list = []
    for base in seq:
        comp_seq_list.append(COMPLEMENTARY_BASE_PAIRS[base])

    comp_seq = ''.join(comp_seq_list)
    return comp_seq[::-1]


# Read sequences from file
scaffold = read_seqfile(SCAFFOLD_SEQFILE)[0]
staples = read_seqfile(STAPLE_SEQFILE)

# Find starting base indices in scaffold
scaffold_domains = []
staple_domains = []
scaffold_identities = []
staple_identities = []
scaffold_indices = []
identity = 1
for staple in staples:
    reverse_comp_staple = reverse_complement(staple)
    if len(staple) in [15, 16]:
        domain_identity = identity
        identity += 1
        scaffold_domain = reverse_comp_staple
        domain_index = scaffold.find(scaffold_domain)
        scaffold_domains.append(scaffold_domain)
        scaffold_indices.append(domain_index)
        scaffold_identities.append(domain_identity)

        staple_domain = reverse_complement(scaffold_domain)
        staple_domains.append([staple_domain])
        staple_identities.append([domain_identity])
        continue

    domain_1_identity = identity
    identity += 1
    domain_2_identity = identity
    identity += 1
    if len(staple) == 32:
        scaffold_domain_1 = reverse_comp_staple[:16]
        domain_1_index = scaffold.find(scaffold_domain_1)
        scaffold_domain_2 = reverse_comp_staple[16:32]
        domain_2_index = scaffold.find(scaffold_domain_2)
    elif len(staple) == 31:
        scaffold_domain_1 = reverse_comp_staple[:16]
        domain_1_index = scaffold.find(scaffold_domain_1)
        if domain_1_index == -1:
            scaffold_domain_1 = reverse_comp_staple[:15]
            domain_1_index = scaffold.find(scaffold_domain_1)
            scaffold_domain_2 = reverse_comp_staple[15:31]
            domain_2_index = scaffold.find(scaffold_domain_2)
        else:
            domain_1_index = scaffold.find(scaffold_domain_1)
            scaffold_domain_2 = reverse_comp_staple[16:31]
            domain_2_index = scaffold.find(scaffold_domain_2)
    else:
        raise SystemExit

    scaffold_domains.extend([scaffold_domain_1, scaffold_domain_2])
    scaffold_indices.extend([domain_1_index, domain_2_index])
    staple_domain_1 = reverse_complement(scaffold_domain_1)
    staple_domain_2 = reverse_complement(scaffold_domain_2)
    staple_domains.append([staple_domain_1, staple_domain_2])
    staple_identities.append([domain_1_identity, domain_2_identity])
    scaffold_identities.extend([domain_1_identity, domain_2_identity])

# Sort scaffold identities based on there order along the scaffold chain
scaffold_identities = [[-identity for (index, identity) in sorted(zip(
    scaffold_indices, scaffold_identities), key=lambda pair: pair[0])]]

# Reorder scaffold domain sequences to be consistent
scaffold_domains_reordered = []
for scaffold_i in scaffold_identities[0]:
    scaffold_i = -scaffold_i - 1
    seq = scaffold_domains[scaffold_i]
    scaffold_domains_reordered.append(seq)

scaffold_domains = scaffold_domains_reordered
sequences = [scaffold_domains] + staple_domains

# Create linear scaffold positions and orientations
positions = []
orientations = []
if CYCLIC:
    position = np.array([0, 0, 0])
    directions = np.array([1, 0, 0]), np.array(
        [0, -1, 0]), np.array([-1, 0, 0]), np.array([0, 1, 0])
    for direction in directions:
        for j in range(len(scaffold_domains) // 4):
            position += direction
            position = position.tolist()
            positions.append(position)
            orientations.append([1, 0, 0])
else:
    position = np.array([-1, 0, 0])
    for j in range(len(scaffold_domains)):
        position += np.array([1, 0, 0])
        position = position.tolist()
        positions.append(position)
        orientations.append([1, 0, 0])

# Output json file
json_origami = {'origami': {}}
json_origami['origami']['cyclic'] = CYCLIC
identities = scaffold_identities + staple_identities
json_origami['origami']['identities'] = identities
json_origami['origami']['sequences'] = sequences
json_origami['origami']['configurations'] = []
json_origami['origami']['configurations'].append({})
json_origami['origami']['configurations'][0]['step'] = 0
json_origami['origami']['configurations'][0]['chains'] = []
json_origami['origami']['configurations'][0]['chains'].append({})
json_origami['origami']['configurations'][0]['chains'][0]['index'] = 0
json_origami['origami']['configurations'][0]['chains'][0]['identity'] = 0
json_origami['origami']['configurations'][0]['chains'][0]['positions'] = positions
json_origami['origami']['configurations'][0]['chains'][0]['orientations'] = orientations

json.dump(json_origami, open(OUTPUT_FILENAME, 'w'),
          indent=4, separators=(',', ': '))
