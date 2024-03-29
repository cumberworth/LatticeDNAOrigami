// ptmc_simulation.cpp

#include <chrono>
#include <filesystem>
#include <sstream>
#include <string>

#include "LatticeDNAOrigami/files.hpp"
#include "LatticeDNAOrigami/ptmc_simulation.hpp"

namespace ptmc {

namespace fs = std::filesystem;

using std::cout;
using std::min;
using std::pair;
using std::string;
using std::chrono::steady_clock;

using files::OrigamiTrajInputFile;
using origami::Chains;

PTGCMCSimulation::PTGCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        GCMCSimulation(origami_system, ops, biases, params),
        m_num_reps {params.m_num_reps},
        m_swaps {params.m_swaps},
        m_max_pt_dur {params.m_max_pt_dur},
        m_exchange_interval {params.m_exchange_interval},
        m_config_output_freq {params.m_configs_output_freq} {

    string string_rank {std::to_string(m_rank)};

    // Update starting configs if restarting
    if (m_params.m_restart_from_config == true) {
        string filename {
                params.m_restart_traj_filebase + "-" + string_rank +
                params.m_restart_traj_postfix};
        OrigamiTrajInputFile traj_inp {filename};
        Chains restart_config {traj_inp.read_config(params.m_restart_step)};
        m_origami_system.set_config(restart_config);
    }

    string output_filebase {params.m_output_filebase + "-" + string_rank};
    m_output_files = simulation::setup_output_files(
            params,
            output_filebase,
            m_origami_system,
            m_ops,
            m_biases,
            m_random_gens);
    m_logging_stream = new ofstream {output_filebase + ".out"};

    // Initialize quantity index to replica index vectors
    if (m_rank == m_master_rep) {

        // For now just take last line of restart file
        if (m_params.m_restart_from_swap == true) {
            string q_to_repi_str;
            string filename {params.m_restart_swap_file};
            fs::path f {filename};
            if (!fs::exists(f)) {
                throw utility::FileError {
                        "Restart swap file " + filename + " does not exist"};
            }
            std::ifstream swap_file {filename};
            string line;

            // Should wrap this in try catch (or make its own file class)
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
}

void PTGCMCSimulation::initialize_swap_file(InputParameters& params) {
    if (m_rank == m_master_rep) {
        m_swapfile.open(params.m_output_filebase + ".swp");
        for (int rep {0}; rep != m_num_reps; rep++) {
            for (auto q_i: m_exchange_q_is) {
                m_swapfile << m_control_qs[q_i][rep];
                m_swapfile << "/";
            }
            m_swapfile << " ";
        }
        m_swapfile << "\n";
    }
}

void PTGCMCSimulation::run() {
    long long int step {0};

    auto start = steady_clock::now();
    for (int swap_i {1}; swap_i != m_swaps + 1; swap_i++) {

        // Update origami system with current replica's quantities
        update_control_qs();

        // Run the simulation
        simulate(m_exchange_interval, step, false, start);
        update_dependent_qs();
        step += m_exchange_interval;

        if (m_rank == m_master_rep) {
            std::chrono::duration<double> dt {(steady_clock::now() - start)};
            if (dt.count() > m_max_pt_dur) {
                master_send_kill(swap_i);
                cout << "Maximum time allowed reached\n";
                break;
            }
        }

        // Send information from slave nodes to master nodes
        if (m_rank != m_master_rep) {
            slave_send(swap_i);
            if (not slave_receive(swap_i)) {
                break;
            }
        }

        // Attempt exchanges between replicas
        else {
            write_swap_entry(step);
            attempt_exchange(swap_i);
        }
    }

    // Write end-of-simulation data
    if (m_rank == m_master_rep) {
        write_swap_entry(step);
        write_acceptance_freqs();
        m_swapfile.close();
    }
}

void PTGCMCSimulation::update_dependent_qs() {
    m_origami_system.update_enthalpy_and_entropy();
    double DH {m_origami_system.hybridization_enthalpy()};
    double D_stacking {m_origami_system.stacking_energy()};
    double bias_e {m_biases.get_total_bias()};

    m_replica_dependent_qs[m_enthalpy_i] = DH;
    m_replica_dependent_qs[m_bias_i] = bias_e;
    m_replica_dependent_qs[m_stacking_i] = D_stacking;
}

void PTGCMCSimulation::slave_send(int swap_i) {
    for (size_t q_i {0}; q_i != m_replica_dependent_qs.size(); q_i++) {
        double q {m_replica_dependent_qs[q_i]};
        m_world.send(m_master_rep, swap_i, q);
    }
    for (auto staple_u: m_origami_system.m_staple_us) {
        m_world.send(m_master_rep, swap_i, staple_u);
    }
    for (auto staple_n: m_origami_system.get_staple_counts()) {
        double staple_n_d {static_cast<double>(staple_n)};
        m_world.send(m_master_rep, swap_i, staple_n_d);
    }
}

bool PTGCMCSimulation::slave_receive(int swap_i) {
    for (auto i: m_exchange_q_is) {
        m_world.recv(m_master_rep, swap_i, m_replica_control_qs[i]);
        if (m_replica_control_qs[i] == 999.0) {
            return false;
        }
    }

    return true;
}

void PTGCMCSimulation::master_receive(
        int swap_i,
        vector<vector<double>>& dependent_qs,
        vector<vector<vector<double>>>& per_staple_dependent_qs) {
    master_get_dependent_qs(dependent_qs, per_staple_dependent_qs);
    size_t num_staple_types {m_origami_system.m_identities.size() - 1};
    for (int rep_i {1}; rep_i != m_num_reps; rep_i++) {
        for (size_t i {0}; i != dependent_qs.size(); i++) {
            double q;
            m_world.recv(rep_i, swap_i, q);
            dependent_qs[i].push_back(q);
        }
        for (size_t i {0}; i != per_staple_dependent_qs.size(); i++) {
            vector<double> q_per_staple {};
            for (size_t j {0}; j != num_staple_types; j++) {
                double q;
                m_world.recv(rep_i, swap_i, q);
                q_per_staple.push_back(q);
            }
            per_staple_dependent_qs[i].push_back(q_per_staple);
        }
    }
}

void PTGCMCSimulation::master_send(int swap_i) {
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

void PTGCMCSimulation::master_send_kill(int swap_i) {
    for (int q_i {0}; q_i != m_num_reps; q_i++) {
        int rep_i {m_q_to_repi[q_i]};
        if (rep_i == m_master_rep) {
            continue;
        }
        else {
            for (size_t i {0}; i != m_exchange_q_is.size(); i++) {
                double msg {999.0};
                m_world.send(rep_i, swap_i, msg);
            }
        }
    }
}

void PTGCMCSimulation::master_get_dependent_qs(
        vector<vector<double>>& dependent_qs,
        vector<vector<vector<double>>>& per_staple_dependent_qs) {
    dependent_qs[m_enthalpy_i].push_back(m_replica_dependent_qs[m_enthalpy_i]);
    dependent_qs[m_bias_i].push_back(m_replica_dependent_qs[m_bias_i]);
    dependent_qs[m_stacking_i].push_back(m_replica_dependent_qs[m_stacking_i]);
    per_staple_dependent_qs[0].push_back(m_origami_system.m_staple_us);
    vector<int> staple_ns {m_origami_system.get_staple_counts()};
    vector<double> staple_ns_d {staple_ns.begin(), staple_ns.end()};
    per_staple_dependent_qs[1].push_back(staple_ns_d);
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

// This does not do the general bias correctly
// I need the bias of both configurations in both states
double PTGCMCSimulation::calc_acceptance_p(
        vector<pair<double, double>> control_q_pairs,
        vector<pair<double, double>> dependent_q_pairs,
        vector<pair<vector<double>, vector<double>>>
                per_staple_dependent_q_pairs) {

    double temp1 {control_q_pairs[m_temp_i].first};
    double temp2 {control_q_pairs[m_temp_i].second};
    double stacking_mult1 {control_q_pairs[m_stacking_mult_i].first};
    double stacking_mult2 {control_q_pairs[m_stacking_mult_i].second};

    double enthalpy1 {dependent_q_pairs[m_enthalpy_i].first};
    double enthalpy2 {dependent_q_pairs[m_enthalpy_i].second};
    double bias1 {dependent_q_pairs[m_bias_i].first};
    double bias2 {dependent_q_pairs[m_bias_i].second};
    double stacking1 {dependent_q_pairs[m_stacking_i].first};
    double stacking2 {dependent_q_pairs[m_stacking_i].second};

    size_t num_staple_types {m_origami_system.m_identities.size()};
    double DBU_DN {0};
    for (size_t i {0}; i != num_staple_types; i++) {
        double N1 {per_staple_dependent_q_pairs[1].first[i]};
        double N2 {per_staple_dependent_q_pairs[1].second[i]};
        double staple_u1 {per_staple_dependent_q_pairs[0].first[i]};
        double staple_u2 {per_staple_dependent_q_pairs[0].second[i]};
        DBU_DN += (staple_u2 / temp2 - staple_u1 / temp1) * (N2 - N1);
    }

    // Energies are actually E/B, so multiply by T
    double DB {1 / temp2 - 1 / temp1};
    double DH {enthalpy2 * temp2 - enthalpy1 * temp1};
    double Dstacking {stacking2 * temp2 - stacking1 * temp1};
    double DBM {stacking_mult2 / temp2 - stacking_mult1 / temp1};
    double DBias {bias2 * temp2 - bias1 * temp1};
    double p_accept {
            min({1.0, exp(DB * (DH + DBias) + DBM * Dstacking - DBU_DN)})};

    return p_accept;
}

void PTGCMCSimulation::write_swap_entry(long long int step) {
    if (step % m_config_output_freq == 0) {
        for (auto repi: m_q_to_repi) {
            m_swapfile << repi << " ";
        }
        m_swapfile << "\n";
    }
}

OneDPTGCMCSimulation::OneDPTGCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        PTGCMCSimulation(origami_system, ops, biases, params),
        m_attempt_count(params.m_temps.size() - 1, 0),
        m_swap_count(params.m_temps.size() - 1, 0) {

    initialize_control_qs(params);
}

// Could probably break this into two methods
void OneDPTGCMCSimulation::initialize_control_qs(InputParameters& params) {

    // Initialize master variables
    if (m_rank == m_master_rep) {
        m_control_qs.push_back(params.m_temps);
        m_control_qs.push_back(params.m_chem_pot_mults);
        m_control_qs.push_back(params.m_bias_mults);
        m_control_qs.push_back(params.m_stacking_mults);
    }

    // Initialize quantities of each replica (updating origami happens in run)
    for (int i {0}; i != m_num_reps; i++) {
        if (m_rank == i) {
            m_replica_control_qs[m_temp_i] = params.m_temps[i];
            m_replica_control_qs[m_staple_u_mult_i] =
                    params.m_chem_pot_mults[i];
            m_replica_control_qs[m_bias_i] = params.m_bias_mults[i];
            m_replica_control_qs[m_stacking_mult_i] =
                    params.m_stacking_mults[i];
        }
    }
}

void OneDPTGCMCSimulation::attempt_exchange(int swap_i) {

    // Collect results from all replicas
    vector<vector<double>> dependent_qs {{}, {}, {}, {}};
    vector<vector<vector<double>>> per_staple_dependent_qs {{}, {}};
    master_receive(swap_i, dependent_qs, per_staple_dependent_qs);

    // Iterate through pairs in current set and attempt swap
    int swap_set {swap_i % 2};
    for (int i {swap_set}; i < (m_num_reps - 1); i += 2) {
        m_attempt_count[i]++;

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
        vector<pair<vector<double>, vector<double>>>
                per_staple_dependent_q_pairs {};
        for (auto per_staple_dependent_q: per_staple_dependent_qs) {
            vector<double> q_1 {per_staple_dependent_q[repi1]};
            vector<double> q_2 {per_staple_dependent_q[repi2]};
            per_staple_dependent_q_pairs.push_back({q_1, q_2});
        }

        double p_accept {calc_acceptance_p(
                control_q_pairs,
                dependent_q_pairs,
                per_staple_dependent_q_pairs)};
        bool accept {test_acceptance(p_accept)};

        // If accepted swap indices from quantities to replicas
        if (accept) {
            m_swap_count[i]++;
            m_q_to_repi[i] = repi2;
            m_q_to_repi[i + 1] = repi1;
        }
    }

    // Send updated temperatures and chem. pots to slaves
    master_send(swap_i);
}

void OneDPTGCMCSimulation::write_acceptance_freqs() {

    for (size_t i {0}; i != m_attempt_count.size(); i++) {
        cout << m_control_qs[m_temp_i][i] << " ";
        cout << m_control_qs[m_temp_i][i + 1] << " ";
        cout << m_swap_count[i] << " ";
        cout << m_attempt_count[i] << " ";
        cout << (static_cast<double>(m_swap_count[i]) / m_attempt_count[i])
             << " ";
        cout << "\n";
    }
    cout << "\n";
}

TwoDPTGCMCSimulation::TwoDPTGCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        PTGCMCSimulation(origami_system, ops, biases, params),
        m_v1_dim {static_cast<int>(params.m_temps.size())},
        m_v2_dim {static_cast<int>(params.m_stacking_mults.size())},
        m_v1s {params.m_temps},
        m_v2s {params.m_stacking_mults},
        m_attempt_count(
                2,
                vector<vector<int>>(m_v1_dim, vector<int>(m_v2_dim, 0))),
        m_swap_count(
                2,
                vector<vector<int>>(m_v1_dim, vector<int>(m_v2_dim, 0))) {

    initialize_control_qs(params);
    m_exchange_q_is.push_back(m_temp_i);
    m_exchange_q_is.push_back(m_staple_u_mult_i);
    m_exchange_q_is.push_back(m_stacking_mult_i);
    initialize_swap_file(params);
}

// Could probably break this into two methods
// This is a hack version that only works with T as v1 and stacking as v2
void TwoDPTGCMCSimulation::initialize_control_qs(InputParameters& params) {

    vector<double> temps {};
    vector<double> staple_us {};
    vector<double> bias_mults {};
    vector<double> stacking_mults {};
    for (size_t v1_i {0}; v1_i != params.m_temps.size(); v1_i++) {
        for (size_t v2_i {0}; v2_i != params.m_stacking_mults.size(); v2_i++) {
            double temp {params.m_temps[v1_i]};
            temps.push_back(temp);
            double staple_u_mult;
            staple_us.push_back(params.m_staple_u_mult);
            bias_mults.push_back(1);
            double stacking_mult {params.m_stacking_mults[v2_i]};
            stacking_mults.push_back(stacking_mult);
        }
    }

    if (m_rank == m_master_rep) {
        //  Temps
        m_control_qs.push_back(temps);

        // Chemical potentials and volumes
        m_control_qs.push_back(staple_us);

        // Biases
        vector<double> bias_mults(temps.size(), 1);
        m_control_qs.push_back(bias_mults);

        // Stacks
        m_control_qs.push_back(stacking_mults);
    }

    // Initialize quantities of each replica (updating on origami happens
    // in run)
    m_replica_control_qs[m_temp_i] = temps[m_rank];
    m_replica_control_qs[m_staple_u_mult_i] = staple_us[m_rank];
    m_replica_control_qs[m_bias_i] = bias_mults[m_rank];
    m_replica_control_qs[m_stacking_mult_i] = stacking_mults[m_rank];
}

void TwoDPTGCMCSimulation::attempt_exchange(int swap_i) {

    // Collect results from all replicas
    vector<vector<double>> dependent_qs {{}, {}, {}, {}};
    vector<vector<vector<double>>> per_staple_dependent_qs {{}, {}};
    master_receive(swap_i, dependent_qs, per_staple_dependent_qs);

    // Iterate through pairs in current set and attempt swap
    int swap_set {swap_i % 4};
    int swap_v {swap_i % 2};
    int i_start {m_i_starts[swap_set]};
    int j_start {m_j_starts[swap_set]};
    int i_incr {m_i_incrs[swap_set]};
    int j_incr {m_j_incrs[swap_set]};
    int i_end {m_i_ends[swap_set]};
    int j_end {m_j_ends[swap_set]};
    int rep_incr {m_rep_incrs[swap_set]};
    for (int i {i_start}; i < (i_end); i += i_incr) {
        for (int j {j_start}; j < (j_end); j += j_incr) {
            int rep_i {i * m_v2_dim + j};
            int rep_j {rep_i + rep_incr};
            m_attempt_count[swap_v][i][j]++;

            // Collect values
            vector<pair<double, double>> control_q_pairs {};
            for (auto control_q: m_control_qs) {
                double q_1 {control_q[rep_i]};
                double q_2 {control_q[rep_j]};
                control_q_pairs.push_back({q_1, q_2});
            }

            int repi1 {m_q_to_repi[rep_i]};
            int repi2 {m_q_to_repi[rep_j]};
            vector<pair<double, double>> dependent_q_pairs {};
            for (auto dependent_q: dependent_qs) {
                double q_1 {dependent_q[repi1]};
                double q_2 {dependent_q[repi2]};
                dependent_q_pairs.push_back({q_1, q_2});
            }
            vector<pair<vector<double>, vector<double>>>
                    per_staple_dependent_q_pairs {};
            for (auto per_staple_dependent_q: per_staple_dependent_qs) {
                vector<double> q_1 {per_staple_dependent_q[repi1]};
                vector<double> q_2 {per_staple_dependent_q[repi2]};
                per_staple_dependent_q_pairs.push_back({q_1, q_2});
            }

            double p_accept {calc_acceptance_p(
                    control_q_pairs,
                    dependent_q_pairs,
                    per_staple_dependent_q_pairs)};

            bool accept {test_acceptance(p_accept)};

            // If accepted swap indices from quantities to replicas
            if (accept) {
                m_swap_count[swap_v][i][j]++;
                m_q_to_repi[rep_i] = repi2;
                m_q_to_repi[rep_j] = repi1;
            }
        }
    }

    // Send updated temperatures and chem. pots to slaves
    master_send(swap_i);
}

void TwoDPTGCMCSimulation::write_acceptance_freqs() {

    for (int v1_i {0}; v1_i != (m_v1_dim - 1); v1_i++) {
        for (int v2_i {0}; v2_i != m_v2_dim; v2_i++) {
            cout << m_v1s[v1_i] << " ";
            cout << m_v1s[v1_i + 1] << " ";
            cout << m_v2s[v2_i] << " ";
            int swap_count {m_swap_count[0][v1_i][v2_i]};
            int attempt_count {m_attempt_count[0][v1_i][v2_i]};
            cout << swap_count << " ";
            cout << attempt_count << " ";
            cout << (static_cast<double>(swap_count) / attempt_count) << " ";
            cout << "\n";
        }
    }
    cout << "\n";

    for (int v2_i {0}; v2_i != (m_v2_dim - 1); v2_i++) {
        for (int v1_i {0}; v1_i != m_v1_dim; v1_i++) {
            cout << m_v2s[v2_i] << " ";
            cout << m_v2s[v2_i + 1] << " ";
            cout << m_v1s[v1_i] << " ";
            int swap_count {m_swap_count[1][v1_i][v2_i]};
            int attempt_count {m_attempt_count[1][v1_i][v2_i]};
            cout << swap_count << " ";
            cout << attempt_count << " ";
            cout << (static_cast<double>(swap_count) / attempt_count) << " ";
            cout << "\n";
        }
    }
    cout << "\n";
}

void TwoDPTGCMCSimulation::update_control_qs() {
    double temp {m_replica_control_qs[m_temp_i]};
    double stacking_mult {m_replica_control_qs[m_stacking_mult_i]};
    m_origami_system.update_temp(temp, stacking_mult);
    double staple_u_mult {m_replica_control_qs[m_staple_u_mult_i]};
    m_origami_system.update_staple_us(temp, staple_u_mult);
}

TPTGCMCSimulation::TPTGCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        OneDPTGCMCSimulation(origami_system, ops, biases, params) {

    m_exchange_q_is.push_back(m_temp_i);
    initialize_swap_file(params);
}

STPTGCMCSimulation::STPTGCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        OneDPTGCMCSimulation(origami_system, ops, biases, params) {

    m_exchange_q_is.push_back(m_temp_i);
    m_exchange_q_is.push_back(m_stacking_mult_i);
    initialize_swap_file(params);
}

UTPTGCMCSimulation::UTPTGCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        OneDPTGCMCSimulation(origami_system, ops, biases, params) {

    m_exchange_q_is.push_back(m_temp_i);
    m_exchange_q_is.push_back(m_staple_u_mult_i);
    initialize_swap_file(params);
}

HUTPTGCMCSimulation::HUTPTGCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        OneDPTGCMCSimulation(origami_system, ops, biases, params) {

    m_exchange_q_is.push_back(m_temp_i);
    m_exchange_q_is.push_back(m_staple_u_mult_i);
    m_exchange_q_is.push_back(m_bias_mult_i);
    initialize_swap_file(params);
}

void TPTGCMCSimulation::update_control_qs() {
    double temp {m_replica_control_qs[m_temp_i]};
    m_origami_system.update_temp(temp);
    double staple_u_mult {m_replica_control_qs[m_staple_u_mult_i]};
    m_origami_system.update_staple_us(temp, staple_u_mult);
}

void UTPTGCMCSimulation::update_control_qs() {
    double temp {m_replica_control_qs[m_temp_i]};
    m_origami_system.update_temp(temp);
    double staple_u_mult {m_replica_control_qs[m_staple_u_mult_i]};
    m_origami_system.update_staple_us(temp, staple_u_mult);
}

void HUTPTGCMCSimulation::update_control_qs() {
    double temp {m_replica_control_qs[m_temp_i]};
    m_origami_system.update_temp(temp);
    double staple_u_mult {m_replica_control_qs[m_staple_u_mult_i]};
    m_origami_system.update_staple_us(temp, staple_u_mult);
    double bias_mult {m_replica_control_qs[m_bias_mult_i]};
    m_origami_system.update_bias_mult(bias_mult);
}

void STPTGCMCSimulation::update_control_qs() {
    double temp {m_replica_control_qs[m_temp_i]};
    double stacking_mult {m_replica_control_qs[m_stacking_mult_i]};
    m_origami_system.update_temp(temp, stacking_mult);
    double staple_u_mult {m_replica_control_qs[m_staple_u_mult_i]};
    m_origami_system.update_staple_us(temp, staple_u_mult);
}
} // namespace ptmc
