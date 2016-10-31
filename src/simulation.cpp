// simulation.cpp

#include <random>
#include <iostream>
#include <set>
#include <algorithm>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "utility.h"
#include "random_gens.h"
#include "movetypes.h"
#include "simulation.h"

using std::cout;
using std::min;

namespace mpi = boost::mpi;

using namespace Movetypes;
using namespace Simulation;
using namespace Utility;
using namespace RandomGen;

GCMCSimulation::GCMCSimulation(OrigamiSystem& origami_system,
        InputParameters params) :
        m_origami_system {origami_system} {

    m_logging_freq = params.m_logging_freq;
    m_centering_freq = params.m_centering_freq;
    m_movetype_constructors = params.m_movetype_constructors;

    // Create cumulative probability array
    double cum_prob {0};
    for (size_t i {0}; i != m_movetype_constructors.size(); i++) {
        cum_prob += params.m_movetype_probs[i];
        m_cumulative_probs.push_back(cum_prob);
    }
}

GCMCSimulation::~GCMCSimulation() {
    for (auto output_file: m_output_files) {
        delete output_file;
    }
}

void GCMCSimulation::simulate(int steps, int start_step) {

    for (int step {start_step + 1}; step != (steps + start_step + 1); step ++) {
        unique_ptr<MCMovetype> movetype {select_movetype()};
        bool accepted;
        accepted = movetype->attempt_move();

        if (not accepted) {
            movetype->reset_origami();
        }

        if (m_centering_freq != 0 and step % m_centering_freq == 0) {
            m_origami_system.centre();
            m_origami_system.check_all_constraints();
        }

        // Write log entry to standard out
        if (m_logging_freq !=0 and step % m_logging_freq == 0) {
            write_log_entry(step, *movetype, accepted);
        }

        // Write system properties to file
        for (auto output_file: m_output_files) {
            if (output_file->m_write_freq != 0 and step % output_file->m_write_freq == 0) {
                output_file->write(step);
            }
        }
    }
}

unique_ptr<MCMovetype> GCMCSimulation::select_movetype() {
    unique_ptr<MCMovetype> movetype;
    double prob {m_random_gens.uniform_real()};
    for (size_t i {0}; i != m_cumulative_probs.size(); i++) {
        if (prob < m_cumulative_probs[i]) {
            movetype = m_movetype_constructors[i](m_origami_system,
                    m_random_gens, m_ideal_random_walks);
            break;
        }
    }
    return movetype;
}

void GCMCSimulation::write_log_entry(int step, MCMovetype& movetype,
        bool accepted) {

    *m_logging_stream << "Step: " << step << " ";
    *m_logging_stream << "Movetype: " << movetype.m_label() << " ";
    *m_logging_stream << "Staples: " << m_origami_system.num_staples() << " ";
    *m_logging_stream << "Accepted: " << accepted << " ";
    *m_logging_stream << "Temp: " << m_origami_system.m_temp << " ";
    *m_logging_stream << "Energy: " << m_origami_system.energy() << " ";
    *m_logging_stream << "\n";
}

ConstantTGCMCSimulation::ConstantTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters params) :
        GCMCSimulation(origami_system, params),
        m_steps {params.m_steps} {

    m_logging_stream = &cout;
    m_output_files = setup_output_files(params, params.m_output_filebase,
            m_origami_system);
}

AnnealingGCMCSimulation::AnnealingGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters params) :
        GCMCSimulation(origami_system, params),
        m_max_temp {params.m_max_temp},
        m_min_temp {params.m_min_temp},
        m_temp_interval {params.m_temp_interval},
        m_steps_per_temp {params.m_steps_per_temp} {

    if (fmod(m_max_temp - m_min_temp, m_temp_interval) != 0) {
        cout << "Bad temperature interval";
    }
    m_logging_stream = &cout;
    m_output_files = setup_output_files(params, params.m_output_filebase,
            m_origami_system);
}

void AnnealingGCMCSimulation::run() {
    double temp {m_max_temp};
    int step {0};
    while (temp >= m_min_temp) {
        m_origami_system.update_temp(temp);
        simulate(m_steps_per_temp, step);
        temp -= m_temp_interval;
        step += m_steps_per_temp;
    }
}

PTGCMCSimulation::PTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters params) :
        GCMCSimulation(origami_system, params),
        m_num_reps {params.m_num_reps},
        m_exchange_interval {params.m_exchange_interval} {

    m_swaps = params.m_steps / m_exchange_interval;
    string string_rank {std::to_string(m_rank)};
    string output_filebase {params.m_output_filebase + "-" + string_rank};
    m_output_files = setup_output_files(params, output_filebase,
            m_origami_system);
    m_logging_stream = new ofstream {output_filebase + ".out"};

    // Initialize temp of each replica
    for (int i {0}; i != m_num_reps; i++) {
        if (m_rank == i) {
            m_temp = params.m_temps[i];
        }
    }

    // Initialize temperature and temperature index to replica index vectors
    if (m_rank == m_master_rep) {
        m_temps = params.m_temps;
        for (int i {0}; i != m_num_reps; i++) {
            m_tempi_to_repi.push_back(i);
        }
    }

    // Initalize swap file
    if (m_rank == m_master_rep) {
        m_swapfile.open(params.m_output_filebase + ".swp");
        for (auto temp: m_temps) {
            m_swapfile <<  temp << " ";
        }
        m_swapfile << "\n";
    }
}

void PTGCMCSimulation::run() {
    int step {0};
    vector<int> attempt_count(m_num_reps - 1, 0);
    vector<int> swap_count(m_num_reps - 1, 0);
    for (int swap_i {1}; swap_i != m_swaps + 1; swap_i++) {
        m_origami_system.update_temp(m_temp);
        simulate(m_exchange_interval, step);
        step += m_exchange_interval;
        if (m_rank != m_master_rep) {
            send_and_recieve_exchange_info(swap_i);
        }
        else {
            attempt_exchange(swap_i, attempt_count, swap_count);
            m_temp = m_temps[m_master_rep];
            write_swap_entry();
        }
    }
    if (m_rank == m_master_rep) {
        write_acceptance_freqs(attempt_count, swap_count);
        m_swapfile.close();
    }
}

void PTGCMCSimulation::send_and_recieve_exchange_info(int swap_i) {

    double energy {m_origami_system.energy()};
    m_world.send(m_master_rep, swap_i, energy);
    m_world.recv(m_master_rep, swap_i, m_temp);
}

void PTGCMCSimulation::attempt_exchange(int swap_i,
        vector<int>& attempt_count, vector<int>& swap_count) {

    // Collect results and add own
    vector<double> energies {m_origami_system.energy()};
    for (int rep_i {1}; rep_i != m_num_reps; rep_i++) {
        double energy;
        m_world.recv(rep_i, swap_i, energy);
        energies.push_back(energy);
    }

    // Iterate through pairs in current set and attempt swap
    int swap_set {swap_i % 2};
    for (int i {swap_set}; i < (m_num_reps - 1); i += 2) {
        attempt_count[i]++;

        // Collect values
        double temp1 {m_temps[i]};
        double temp2 {m_temps[i + 1]};
        int repi1 {m_tempi_to_repi[i]};
        int repi2 {m_tempi_to_repi[i + 1]};

        // Energies are actually E/B, so multiply by T
        double energy1 {energies[repi1] * temp1};
        double energy2 {energies[repi2] * temp2};
        bool accept {test_acceptance(temp1, temp2, energy1, energy2)};
        if (accept) {
            swap_count[i]++;
            int repi1_old {repi1};
            repi1 = repi2;
            repi2 = repi1_old;
            m_tempi_to_repi[i] = repi1;
            m_tempi_to_repi[i + 1] = repi2;
        }
    }

    // Send temps
    for (int rep_i {1}; rep_i != m_num_reps; rep_i++) {
        m_world.send(rep_i, swap_i, m_temps[rep_i]);
    }
}

bool PTGCMCSimulation::test_acceptance(double temp1, double temp2, double energy1,
        double energy2) {

    double DB {1/temp2 - 1/temp1};
    double DE {energy2 - energy1};
    double p_accept {min({1.0, exp(-DB*DE)})};
    bool accept;
    if (p_accept == 1) {
        accept = true;
    }
    else {
        double prob {m_random_gens.uniform_real()};
        if (p_accept > prob) {
            accept = true;
        }
        else {
            accept = false;
        }
    }

    return accept;
}

void PTGCMCSimulation::write_swap_entry() {
    for (auto repi: m_tempi_to_repi) {
        m_swapfile << repi << " ";
    }
    m_swapfile << "\n";
}

void PTGCMCSimulation::write_acceptance_freqs(vector<int> attempt_count,
        vector<int> swap_count) {

    for (size_t i {0}; i != attempt_count.size(); i++) {
        cout << m_temps[i] << " ";
        cout << m_temps[i + 1] << " ";
        cout << swap_count[i] << " ";
        cout << attempt_count[i] << " ";
        cout << (static_cast<double>(swap_count[i]) / attempt_count[i]) << " ";
        cout << "\n";
    }
    cout << "\n";
}

vector<OrigamiOutputFile*> Simulation::setup_output_files(
        InputParameters params, string output_filebase,
        OrigamiSystem& origami) {

    vector<OrigamiOutputFile*> outs {};
    if (params.m_configs_output_freq != 0) {
        OrigamiOutputFile* config_out = new OrigamiTrajOutputFile {
                output_filebase + ".trj", params.m_configs_output_freq,
                origami};
        outs.push_back(config_out);
        OrigamiOutputFile* counts_out = new OrigamiCountsOutputFile {
                output_filebase + ".counts", params.m_counts_output_freq,
                origami};
        outs.push_back(counts_out);
    }

    return outs;
}
