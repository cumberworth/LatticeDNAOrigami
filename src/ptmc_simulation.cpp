// ptmc_simulation.cpp

#include <sstream>

#include "ptmc_simulation.h"

using std::min;

using namespace Simulation;
using namespace PTMC;

PTGCMCSimulation::PTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        GCMCSimulation(origami_system, params),
        m_num_reps {params.m_num_reps},
        m_exchange_interval {params.m_exchange_interval} {
    m_swaps = params.m_steps / m_exchange_interval;
    string string_rank {std::to_string(m_rank)};

    // Update starting configs if restarting
    if (m_params.m_restart_traj_filebase != "") {
        string filename {params.m_restart_traj_filebase + "-" + string_rank +
                ".trj"};
        OrigamiTrajInputFile traj_inp {filename};
        Chains restart_config {traj_inp.read_config(params.m_restart_step)};
        m_origami_system.set_config(restart_config);
    }

    string output_filebase {params.m_output_filebase + "-" + string_rank};
    m_output_files = setup_output_files(params, output_filebase,
            m_origami_system);
    m_logging_stream = new ofstream {output_filebase + ".out"};

    // Initialize quantity index to replica index vectors
    if (m_rank == m_master_rep) {

        // For now just take last line of restart file
        if (params.m_restart_swap_file != "") {
            string q_to_repi_str;
            std::ifstream swap_file {params.m_restart_swap_file};
            string line;
            while (std::getline(swap_file, line)) {
                q_to_repi_str = line;
            }
            std::stringstream q_to_repi_sstream {q_to_repi_str};
            for (int i {0}; i != m_num_reps; i++) {
                int repi;
                q_to_repi_sstream >> repi;
                m_q_to_repi.push_back(repi);
            }
        }
        else {
            for (int i {0}; i != m_num_reps; i++) {
                m_q_to_repi.push_back(i);
            }
        }
    }

    initialize_control_qs(params);
}

void PTGCMCSimulation::initialize_swap_file(InputParameters& params) {
    if (m_rank == m_master_rep) {
        m_swapfile.open(params.m_output_filebase + ".swp");
        for (int rep {0}; rep != m_num_reps; rep++) {
            for (auto q_i: m_exchange_q_is) {
                m_swapfile <<  m_control_qs[q_i][rep];
                if (q_i + 1 != static_cast<int>(m_exchange_q_is.size())) {
                    m_swapfile << "/";
                }
                else {
                    m_swapfile << " ";
                }
            }
        }
        m_swapfile << "\n";
        write_swap_entry();
    }
}

PTGCMCSimulation::~PTGCMCSimulation() {
    delete m_logging_stream;
}

// Could probably break this into two methods
void PTGCMCSimulation::initialize_control_qs(InputParameters& params) {
    if (m_rank == m_master_rep) {
 
        //  Temps
        m_control_qs.push_back(params.m_temps);

        // Chemical potentials
        m_control_qs.push_back({});
        for (int i {0}; i != m_num_reps; i++) {

            // Calculate chemical potential of each replica if constant [staple]
            double staple_u;
            if (m_params.m_constant_staple_M) {
                staple_u = molarity_to_chempot(m_params.m_staple_M, params.m_temps[i],
                        params.m_lattice_site_volume);
            }
            else {
                staple_u = molarity_to_chempot(m_params.m_staple_M,
                        params.m_temp_for_staple_u,
                params.m_lattice_site_volume);
                staple_u *= m_params.m_chem_pot_mults[i];
            }
            m_control_qs[m_staple_u_i].push_back(staple_u);
        }

        // Biases
        m_control_qs.push_back(params.m_bias_mults);
    }

    // Initialize quantities of each replica (updating on origami happens in run)
    for (int i {0}; i != m_num_reps; i++) {
        if (m_rank == i) {
            m_replica_control_qs[m_temp_i] = params.m_temps[i];
            m_replica_control_qs[m_bias_i] = params.m_bias_mults[i];

            // Update chemical potential of each replica if constant [staple]
            // Recalculating for each node rather than sending from master
            if (m_params.m_constant_staple_M) {
                m_replica_control_qs[m_staple_u_i] = molarity_to_chempot(
                        m_params.m_staple_M, m_replica_control_qs[m_temp_i],
                        params.m_lattice_site_volume);
            }
            else {
                double staple_u {molarity_to_chempot(m_params.m_staple_M,
                        params.m_temp_for_staple_u,
                        params.m_lattice_site_volume)};
                staple_u *= m_params.m_chem_pot_mults[i];
                m_replica_control_qs[m_staple_u_i] = staple_u;
            }
        }
    }
}

void PTGCMCSimulation::run() {
    // Run temperature/chemical potential parallel tempering grand cannonical monte carlo
    int step {0};

    // Keep track of number of attempted and succesful swaps
    vector<int> attempt_count(m_num_reps - 1, 0);
    vector<int> swap_count(m_num_reps - 1, 0);
    for (int swap_i {1}; swap_i != m_swaps + 1; swap_i++) {

        // Update origami system with current replica's quantities
        update_control_qs();

        // Run the simulation
        simulate(m_exchange_interval, step);
        step += m_exchange_interval;

        // Send information from slave nodes to master nodes
        if (m_rank != m_master_rep) {
            slave_send(swap_i);
            slave_receive(swap_i);
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

void PTGCMCSimulation::update_dependent_qs() {
    double DH {m_origami_system.enthalpy_and_entropy().enthalpy};
    double N {static_cast<double>(m_origami_system.num_staples())};
    double bias_e {m_origami_system.bias()};

    m_replica_dependent_qs[m_enthalpy_i] = DH;
    m_replica_dependent_qs[m_staples_i] = N;
    m_replica_dependent_qs[m_bias_i] = bias_e;
}

void PTGCMCSimulation::slave_send(int swap_i) {
    // Send quantities to master
    for (size_t q_i {0}; q_i != m_replica_dependent_qs.size(); q_i++) {
        double q {m_replica_dependent_qs[q_i]};
        m_world.send(m_master_rep, swap_i, q);
    }
}

void PTGCMCSimulation::slave_receive(int swap_i) {
    // Receive quantities from master
    for (auto i: m_exchange_q_is) {
        m_world.recv(m_master_rep, swap_i, m_replica_control_qs[i]);
    }
}

void PTGCMCSimulation::master_receive(int swap_i,
        vector<vector<double>>& dependent_qs) {
    master_get_dependent_qs(dependent_qs);
    for (int rep_i {1}; rep_i != m_num_reps; rep_i++) {
        for (size_t i {0}; i != dependent_qs.size(); i++) {
            double q;
            m_world.recv(rep_i, swap_i, q);
            dependent_qs[i].push_back(q);
        }
    }
}

void PTGCMCSimulation::master_send(int swap_i) {
    // Send temperatures and chemical potentials to slaves
    for (int q_i {0}; q_i != m_num_reps; q_i++) {
        int rep_i {m_q_to_repi[q_i]};
        if (rep_i == m_master_rep) {
            for (auto i: m_exchange_q_is) {
                m_replica_control_qs[i] = m_control_qs[i][q_i];
            }
        }
        else {
            for (auto i: m_exchange_q_is) {
                m_world.send(rep_i, swap_i, m_control_qs[i][q_i]);
            }
        }
    }
}

void PTGCMCSimulation::master_get_dependent_qs(
        vector<vector<double>>& dependent_qs) {
    double DH {m_origami_system.enthalpy_and_entropy().enthalpy};
    double staples {static_cast<double>(m_origami_system.num_staples())};
    double bias {m_origami_system.bias()};
    dependent_qs[m_enthalpy_i].push_back(DH);
    dependent_qs[m_staples_i].push_back(staples);
    dependent_qs[m_bias_i].push_back(bias);
}

void PTGCMCSimulation::attempt_exchange(int swap_i,
        vector<int>& attempt_count, vector<int>& swap_count) {
    // Alternates between two sets of pairs, allows for changes in T, u, and N

    // Collect results from all replicas
    vector<vector<double>> dependent_qs {{}, {}, {}};
    master_receive(swap_i, dependent_qs);

    // Iterate through pairs in current set and attempt swap
    int swap_set {swap_i % 2};
    for (int i {swap_set}; i < (m_num_reps - 1); i += 2) {
        attempt_count[i]++;

        // Collect values
        vector<pair<double, double>> control_q_pairs {};
        for (auto control_q: m_control_qs) {
            double q_1 {control_q[i]};
            double q_2 {control_q[i + 1]};
            control_q_pairs.push_back({q_1, q_2});
        }

        int repi1 {m_q_to_repi[i]};
        int repi2 {m_q_to_repi[i + 1]};
        vector<pair<double, double>> dependent_q_pairs {};
        for (auto dependent_q: dependent_qs) {
            double q_1 {dependent_q[repi1]};
            double q_2 {dependent_q[repi2]};
            dependent_q_pairs.push_back({q_1, q_2});
        }

        double p_accept {calc_acceptance_p(control_q_pairs, dependent_q_pairs)};
        bool accept {test_acceptance(p_accept)};

        // If accepted swap indices from quantities to replicas
        if (accept) {
            swap_count[i]++;
            m_q_to_repi[i] = repi2;
            m_q_to_repi[i + 1] = repi1;
        }
    }

    // Send updated temperatures and chem. pots to slaves
    master_send(swap_i);
}

bool PTGCMCSimulation::test_acceptance(double p_accept) {
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

double PTGCMCSimulation::calc_acceptance_p(
        vector<pair<double, double>> control_q_pairs,
        vector<pair<double, double>> dependent_q_pairs) {

    double temp1 {control_q_pairs[m_temp_i].first};
    double temp2 {control_q_pairs[m_temp_i].second};
    double staple_u1 {control_q_pairs[m_staple_u_i].first};
    double staple_u2 {control_q_pairs[m_staple_u_i].second};

    double enthalpy1 {dependent_q_pairs[m_enthalpy_i].first};
    double enthalpy2 {dependent_q_pairs[m_enthalpy_i].second};
    double N1 {dependent_q_pairs[m_staples_i].first};
    double N2 {dependent_q_pairs[m_staples_i].second};
    double bias1 {dependent_q_pairs[m_bias_i].first};
    double bias2 {dependent_q_pairs[m_bias_i].second};

    // Energies are actually E/B, so multiply by T
    double DB {1/temp2 - 1/temp1};
    double DH {enthalpy2*temp2 - enthalpy1*temp1};
    double DBias {bias2*temp2 - bias1*temp1};
    double DN {N2 - N1};
    double DBU {staple_u2 / temp2 - staple_u1 / temp1};
    double p_accept {min({1.0, exp(DB*(DH + DBias) - DBU*DN)})};

    return p_accept;
}

void PTGCMCSimulation::write_swap_entry() {
    for (auto repi: m_q_to_repi) {
        m_swapfile << repi << " ";
    }
    m_swapfile << "\n";
}

void PTGCMCSimulation::write_acceptance_freqs(vector<int> attempt_count,
        vector<int> swap_count) {

    for (size_t i {0}; i != attempt_count.size(); i++) {
        cout << m_control_qs[m_temp_i][i] << " ";
        cout << m_control_qs[m_temp_i][i + 1] << " ";
        cout << swap_count[i] << " ";
        cout << attempt_count[i] << " ";
        cout << (static_cast<double>(swap_count[i]) / attempt_count[i]) << " ";
        cout << "\n";
    }
    cout << "\n";
}

TPTGCMCSimulation::TPTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        PTGCMCSimulation(origami_system, params) {
    m_exchange_q_is.push_back(m_temp_i);
    initialize_swap_file(params);
}

UTPTGCMCSimulation::UTPTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        PTGCMCSimulation(origami_system, params) {
    m_exchange_q_is.push_back(m_temp_i);
    m_exchange_q_is.push_back(m_staple_u_i);
    initialize_swap_file(params);
}

HUTPTGCMCSimulation::HUTPTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        PTGCMCSimulation(origami_system, params) {
    m_exchange_q_is.push_back(m_temp_i);
    m_exchange_q_is.push_back(m_staple_u_i);
    m_exchange_q_is.push_back(m_bias_mult_i);
    initialize_swap_file(params);
}

void TPTGCMCSimulation::update_control_qs() {
    double temp {m_replica_control_qs[m_temp_i]};
    m_origami_system.update_temp(temp);
}

void UTPTGCMCSimulation::update_control_qs() {
    double temp {m_replica_control_qs[m_temp_i]};
    m_origami_system.update_temp(temp);
    double staple_u {m_replica_control_qs[m_staple_u_i]};
    m_origami_system.update_staple_u(staple_u);
}

void HUTPTGCMCSimulation::update_control_qs() {
    double temp {m_replica_control_qs[m_temp_i]};
    m_origami_system.update_temp(temp);
    double staple_u {m_replica_control_qs[m_staple_u_i]};
    m_origami_system.update_staple_u(staple_u);
    double bias_mult {m_replica_control_qs[m_bias_mult_i]};
    m_origami_system.update_bias_mult(bias_mult);
}
