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
        SystemBias& system_bias, InputParameters& params) :
        m_origami_system {origami_system},
        m_system_bias {system_bias},
        m_params {params} {

    m_logging_freq = params.m_logging_freq;
    m_centering_freq = params.m_centering_freq;

    // Create movetype constructors
    for (auto movetype_i: params.m_movetypes) {
        m_movetype_constructors.push_back(movetype[movetype_i]);
    }

    // Create cumulative probability array
    double cum_prob {0};
    for (size_t i {0}; i != m_movetype_constructors.size(); i++) {
        cum_prob += params.m_movetype_probs[i];
        m_cumulative_probs.push_back(cum_prob);
    }

    if (params.m_num_walks_filename.size() != 0) {
        std::ifstream num_walks_file {params.m_num_walks_filename};
        boost::archive::binary_iarchive num_walks_arch {num_walks_file};
        num_walks_arch >> m_ideal_random_walks;
    }
}

GCMCSimulation::~GCMCSimulation() {
    for (auto output_file: m_output_files) {
        delete output_file;
    }
}

void GCMCSimulation::simulate(long int steps, int start_step) {

    for (long int step {start_step + 1}; step != (steps + start_step + 1); step ++) {
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
                    m_system_bias, m_random_gens, m_ideal_random_walks, m_params);
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
    *m_logging_stream << "Bias: " << m_system_bias.calc_bias() << " ";
    *m_logging_stream << "\n";
}

ConstantTGCMCSimulation::ConstantTGCMCSimulation(OrigamiSystem& origami_system,
        SystemBias& system_bias, InputParameters& params) :
        GCMCSimulation(origami_system, system_bias, params),
        m_steps {params.m_steps} {

    m_logging_stream = &cout;
    m_output_files = setup_output_files(params, params.m_output_filebase,
            m_origami_system);
}

AnnealingGCMCSimulation::AnnealingGCMCSimulation(OrigamiSystem& origami_system,
        SystemBias& system_bias, InputParameters& params) :
        GCMCSimulation(origami_system, system_bias, params),
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
        SystemBias& system_bias, InputParameters& params) :
        GCMCSimulation(origami_system, system_bias, params),
        m_num_reps {params.m_num_reps},
        m_exchange_interval {params.m_exchange_interval} {

    m_swaps = params.m_steps / m_exchange_interval;
    string string_rank {std::to_string(m_rank)};
    string output_filebase {params.m_output_filebase + "-" + string_rank};
    m_output_files = setup_output_files(params, output_filebase,
            m_origami_system);
    m_logging_stream = new ofstream {output_filebase + ".out"};

    // Initialize temp, u, bias mult, and temperature index to replica index vectors
    if (m_rank == m_master_rep) {
        m_temps = params.m_temps;
        m_bias_mults = params.m_bias_mults;
        for (int i {0}; i != m_num_reps; i++) {
            m_tempi_to_repi.push_back(i);

            // Calculate chemical potential of each replica if constant [staple]
            double staple_u;
            if (m_params.m_constant_staple_M) {
                staple_u = molarity_to_chempot(m_params.m_staple_M, m_temps[i],
                        params.m_lattice_site_volume);
            }
            else {
                staple_u = molarity_to_chempot(m_params.m_staple_M,
                        params.m_temp_for_staple_u,
                params.m_lattice_site_volume);
                staple_u *= m_params.m_chem_pot_mults[i];
            }
            m_staple_us.push_back(staple_u);
        }
    }

    // Initialize temp, u, and bias mult of each replica (updating on origami happens in run)
    for (int i {0}; i != m_num_reps; i++) {
        if (m_rank == i) {
            m_temp = params.m_temps[i];
            m_bias_mult = params.m_bias_mults[i];

            // Update chemical potential of each replica if constant [staple]
            // Recalculating for each node rathr than sending from master
            if (m_params.m_constant_staple_M) {
                m_staple_u = molarity_to_chempot(m_params.m_staple_M, m_temp,
                        params.m_lattice_site_volume);
            }
            else {
                m_staple_u = molarity_to_chempot(m_params.m_staple_M,
                        params.m_temp_for_staple_u,
                        params.m_lattice_site_volume);
                m_staple_u *= m_params.m_chem_pot_mults[i];
            }
        }
    }


    // Initalize swap file
    if (m_rank == m_master_rep) {
        m_swapfile.open(params.m_output_filebase + ".swp");
        for (auto temp: m_temps) {
            m_swapfile <<  temp << " ";
        }
        m_swapfile << "\n";
        write_swap_entry();
    }
}

void PTGCMCSimulation::run() {
    // Run temperature/chemical potential parallel tempering grand cannonical monte carlo
    int step {0};

    // Keep track of number of attempted and succesful swaps
    vector<int> attempt_count(m_num_reps - 1, 0);
    vector<int> swap_count(m_num_reps - 1, 0);
    for (int swap_i {1}; swap_i != m_swaps + 1; swap_i++) {

        // Update origami system with current replica's temp. and chem. pot.
        m_origami_system.update_temp(m_temp);
        m_origami_system.update_staple_u(m_staple_u);

        // Update system bias with current replica's bias multiplier
        m_system_bias.update_bias_mult(m_bias_mult);

        // Run the simulation
        simulate(m_exchange_interval, step);
        step += m_exchange_interval;

        // Send information from slave nodes to master nodes
        if (m_rank != m_master_rep) {
            slave_send_and_recieve(swap_i);
        }

        // Attempt exchanges between replicas
        else {
            attempt_exchange(swap_i, attempt_count, swap_count);
            write_swap_entry();
        }
    }

    // Write end-of-simulation data
    if (m_rank == m_master_rep) {
        write_acceptance_freqs(attempt_count, swap_count);
        m_swapfile.close();
    }
}

void PTGCMCSimulation::slave_send_and_recieve(int swap_i) {
    // Send enthalpy, bias, num staples to master and receive updated temp and u

    ThermoOfHybrid DH_DS {m_origami_system.enthalpy_and_entropy()};
    double bias {m_system_bias.calc_bias()};
    int N {m_origami_system.num_staples()};

    // Send enthalpy, bias, and number of particles to master rep
    m_world.send(m_master_rep, swap_i, DH_DS.enthalpy);
    m_world.send(m_master_rep, swap_i, bias);
    m_world.send(m_master_rep, swap_i, N);

    // Receive temperature and chemical potential from master rep
    m_world.recv(m_master_rep, swap_i, m_temp);
    m_world.recv(m_master_rep, swap_i, m_bias_mult);
    m_world.recv(m_master_rep, swap_i, m_staple_u);
}

void PTGCMCSimulation::master_receive(int swap_i, vector<double>& enthalpies,
        vector<double>& biases, vector<int>& staples) {
    // Receive enthalpies, bias, and number of staples from slaves
    for (int rep_i {1}; rep_i != m_num_reps; rep_i++) {
        double DH;
        double bias;
        int N;
        m_world.recv(rep_i, swap_i, DH);
        m_world.recv(rep_i, swap_i, bias);
        m_world.recv(rep_i, swap_i, N);

        // Add master values
        enthalpies.push_back(DH);
        biases.push_back(bias);
        staples.push_back(N);
    }
}

void PTGCMCSimulation::master_send(int swap_i) {
    // Send temperatures and chemical potentials to slaves
    for (int temp_i {0}; temp_i != m_num_reps; temp_i++) {
        int rep_i {m_tempi_to_repi[temp_i]};
        if (rep_i == m_master_rep) {
            m_temp = m_temps[temp_i];
            m_bias_mult = m_bias_mults[temp_i];
            m_staple_u = m_staple_us[temp_i];
        }
        else {
            m_world.send(rep_i, swap_i, m_temps[temp_i]);
            m_world.send(rep_i, swap_i, m_bias_mults[temp_i]);
            m_world.send(rep_i, swap_i, m_staple_us[temp_i]);
        }
    }
}

void PTGCMCSimulation::attempt_exchange(int swap_i,
        vector<int>& attempt_count, vector<int>& swap_count) {
    // Alternates between two sets of pairs, allows for changes in T, u, and N

    // Collect results from slaves
    ThermoOfHybrid DH_DS {m_origami_system.enthalpy_and_entropy()};
    vector<double> enthalpies {DH_DS.enthalpy};
    vector<double> biases {m_system_bias.calc_bias()};
    vector<int> staples {m_origami_system.num_staples()};
    master_receive(swap_i, enthalpies, biases, staples);

    // Iterate through pairs in current set and attempt swap
    int swap_set {swap_i % 2};
    for (int i {swap_set}; i < (m_num_reps - 1); i += 2) {
        attempt_count[i]++;

        // Collect values
        double temp1 {m_temps[i]};
        double temp2 {m_temps[i + 1]};
        double staple_u1 {m_staple_us[i]};
        double staple_u2 {m_staple_us[i + 1]};
        int repi1 {m_tempi_to_repi[i]};
        int repi2 {m_tempi_to_repi[i + 1]};

        // Energies are actually E/B, so multiply by T
        double enthalpy1 {enthalpies[repi1] * temp1};
        double enthalpy2 {enthalpies[repi2] * temp2};
        double bias1 {enthalpies[repi1] * temp1};
        double bias2 {enthalpies[repi2] * temp2};
        int N1 {staples[repi1]};
        int N2 {staples[repi2]};
        bool accept {test_acceptance(temp1, temp2, staple_u1, staple_u2,
                enthalpy1, enthalpy2, bias1, bias2, N1, N2)};

        // If accepted swap temperatures and chem. pot.s
        if (accept) {
            swap_count[i]++;
            int repi1_old {repi1};
            repi1 = repi2;
            repi2 = repi1_old;
            m_tempi_to_repi[i] = repi1;
            m_tempi_to_repi[i + 1] = repi2;
        }
    }

    // Send updated temperatures and chem. pots to slaves
    master_send(swap_i);
}

bool PTGCMCSimulation::test_acceptance(double temp1, double temp2,
        double staple_u1, double staple_u2, double enthalpy1, double enthalpy2,
        double bias1, double bias2, int N1, int N2) {

    double DB {1/temp2 - 1/temp1};
    double DH {enthalpy2 - enthalpy1};
    double DBias {bias2 - bias1};
    int DN {N2 - N1};
    double DBU {staple_u2 / temp2 - staple_u1 / temp1};
    double p_accept {min({1.0, exp(DB*(DH + DBias) - DBU*DN)})};
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
        InputParameters& params, string output_filebase,
        OrigamiSystem& origami) {
    // Setup trajectory and staple/domain count file

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
