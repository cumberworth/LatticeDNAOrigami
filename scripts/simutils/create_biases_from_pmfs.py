#!/usr/bin/env python3

"""Create bias input files from calculated pmfs"""

import numpy as np
import json


def main():
    filebase = 'inps/snodin-staple_temp-344_run-0_rep-0'
    wins_filename = 'inps/snodin_344_staple.windows'
    wins = read_windows(wins_filename)
    pmfs_filename = 'outs/snodin-long_temp-344_run-1_rep-0_pmfs.sds'
    pmf_array = np.loadtxt(pmfs_filename, skiprows=1)
    pmfs = {(i[0], i[1]): i[2] for i in pmf_array}
    for win in wins:
        biases = {'biases':[]}
        domain_lims = [win[0][0], win[1][0]]
        staple_lims = [win[0][1], win[1][1]]
        for domain in range(domain_lims[0], domain_lims[1] + 1):
            for staple in range(staple_lims[0], staple_lims[1] + 1):
                point = (staple, domain)
                if point not in pmfs:
                    continue
                elif np.isnan(pmfs[point]):
                    continue
                else:
                    biases_entry = {}
                    biases_entry['point'] = [domain, staple]
                    biases_entry['bias'] = -pmfs[point]
                    biases['biases'].append(biases_entry)

        win_filename = create_win_filename(win, filebase)
        json.dump(biases, open(win_filename, 'w'), indent=4,
                separators=(',', ': '))


if __name__ == '__main__':
    main()
