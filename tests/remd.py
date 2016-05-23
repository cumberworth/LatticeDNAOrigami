#!/usr/bin/env python

"""REMD simulation template."""

import multiprocessing as multi
import sys
import time
#sys.path.insert(0, '../../lib/lattice_origami_domains')
sys.path.insert(0, '../')

from lattice_dna_origami.lattice_origami_domains import *

def replica_sim(swap_freq, config_write_freq, count_write_freq, strand_M,
        cation_M, move_settings, input_filename, step, output_filename, temp,
        pipe, rep):

    # Simulation setup
    input_file = JSONInputFile(input_filename)
    origami_system = OrigamiSystemSixteen(input_file, step, temp, strand_M, cation_M)
    energy_lists = {temp: origami_system._hybridization_energies}
    output_file = HDF5OutputFile(output_filename, origami_system,
            config_write_freq=config_write_freq,
            count_write_freq=count_write_freq)
    sim = GCMCSimulation(origami_system, move_settings, output_file)

    while True:

        # Wait to be given a temperature
        pipe.poll(timeout=None)
        temp = pipe.recv()

        # Update origami_system's energy list
        try:
            hybridization_energies = energy_lists[temp]
        except KeyError:
            hybridization_energies = []
            for sequence in origami_system.sequences:
                energy = calc_hybridization_energy(sequence, temp, cation_M)
                hybridization_energies.append(energy)

        sim._origami_system.temp = temp
        sim._origami_system._hybridization_energyes = hybridization_energies

        # Run simulation
        sim.run(swap_freq, logging=0)

        # Send system energy
        senergy = sim._origami_system.energy
        pipe.send(senergy)

def exchange_accepted(energy1, energy2, temp1, temp2):
    DB = (1/temp2) - (1/temp1)
    DE = energy2 - energy1
    ratio = math.exp(DB * DE)
    p_accept = min(1, ratio)
    if p_accept == 1:
        accept = True
    else:
        if p_accept > random.random():
            accept = True
        else:
            accept = False

    return accept

def replica_exchange(pipes, reps, temps, swaps, replicas):
    pair = 0

    # Indexed by temp index, contains indices to reps
    tempi_to_repi = list(range(reps))
    for swap in range(swaps):

        # Wait for reps to finish
        energies = []
        for pipe in pipes:
            pipe.poll(timeout=None)
            energies.append(pipe.recv())

        # Deal with ends that have no partner
        if reps % 2 != 0:
            if pair == 0:
                repi = tempi_to_repi[0]
                pipes[repi].send(temps[0])
            else:
                repi = tempi_to_repi[-1]
                pipes[repi].send(temps[-1])
        else:
            if pair == 1:
                repi = tempi_to_repi[0]
                pipes[repi].send(temps[0])
                repi = tempi_to_repi[-1]
                pipes[repi].send(temps[-1])
            else:
                pass

        # Attempt exchanges on pairs
        for i in range(pair, reps - pair, 2):
            temp1 = temps[i]
            temp2 = temps[i + 1]
            repi1 = tempi_to_repi[i]
            repi2 = tempi_to_repi[i + 1]
            energy1 = energies[repi1]
            energy2 = energies[repi2]
            if exchange_accepted(energy1, energy2, temp1, temp2):
                repi1, repi2 = repi2, repi1
                tempi_to_repi[i] = repi1
                tempi_to_repi[i + 1] = repi2
            else:
                pass

            pipes[repi1].send(temp1)
            pipes[repi2].send(temp2)

        # Switch pair
        if pair == 0:
            pair = 1
        else:
            pair = 0
            
        # Print replica temperatures
        print(*tempi_to_repi, sep=' ')

    # Kill all replicas (this a very harsh method)
    for replica in replicas:
        replica.terminate()

def main():
    reps = 4
    temps = [320, 330, 340, 350]

    swap_freq = 10000
    prod_steps = 100000
    if not value_is_multiple(prod_steps, swap_freq):
        print('Number of production steps must be a multiple of swap frequency.')
        sys.exit()

    swaps = prod_steps // swap_freq

    config_write_freq = 10000
    count_write_freq = 0
    #if not value_is_multiple(swap_freq, config_write_freq) or not (
    #        value_is_multiple(swap_freq, count_write_freq)):
    #    print('Write swap frequencies must be multiples of write frequencies.')
    #    sys.exit()

    strand_M = 1e-3
    cation_M = 1
    move_settings = {MOVETYPE.EXCHANGE_STAPLE: 0.25,
                     MOVETYPE.CB_REGROW_STAPLE: 0.25,
                     MOVETYPE.CB_CONSERVED_TOPOLOGY: 0.25,
                     MOVETYPE.ROTATE_ORIENTATION_VECTOR: 0.25}

    input_filename = 'snodin_unbound.json'
    step = 0
    output_root = 'remd_test'

    replicas = []
    pipes_masterside = []
    for rep, temp in enumerate(temps):
        pipe_masterside, pipe_replicaside = multi.Pipe()
        pipes_masterside.append(pipe_masterside)
        output_filename = output_root + '_{}'.format(rep)

        replica = multi.Process(target=replica_sim, args=(swap_freq,
            config_write_freq, count_write_freq, strand_M, cation_M,
            move_settings, input_filename, step, output_filename, temp,
            pipe_replicaside, rep))

        replica.start()
        pipe_masterside.send(temp)
        replicas.append(replica)

    print(*temps, sep=' ')
    replica_exchange(pipes_masterside, reps, temps, swaps, replicas)

if __name__ == '__main__':
    main()
