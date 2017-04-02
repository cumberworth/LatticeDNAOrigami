// simulation.cpp

#include <random>
#include <iostream>
#include <set>
#include <algorithm>
#include <numeric>
#include <cmath>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "ceres/ceres.h"
#include "glog/logging.h"

#include "utility.h"
#include "random_gens.h"
#include "movetypes.h"
#include "simulation.h"

using std::cout;
using std::min;
using std::pow;

namespace mpi = boost::mpi;

using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::LossFunction;
using ceres::ScaledLoss;

using namespace Movetypes;
using namespace Simulation;
using namespace Utility;
using namespace RandomGen;

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

GCMCSimulation::GCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        m_origami_system {origami_system},
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

    // Load precalculated ideal random walk count data
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
        
        // Pick movetype and apply
        unique_ptr<MCMovetype> movetype {select_movetype()};
        bool accepted;
        accepted = movetype->attempt_move();
        if (not accepted) {
            movetype->reset_origami();
        }

        // Center and check constraints
        if (m_centering_freq != 0 and step % m_centering_freq == 0) {
            m_origami_system.centre();
            m_origami_system.check_all_constraints();
        }

        // Write log entry to standard out
        if (m_logging_freq !=0 and step % m_logging_freq == 0) {
            write_log_entry(step, *movetype, accepted);
        }

        // Update internal simulation variables
        update_internal();

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
                    m_random_gens, m_ideal_random_walks, m_params);
            break;
        }
    }
    return movetype;
}

void GCMCSimulation::write_log_entry(long int step, MCMovetype& movetype,
        bool accepted) {

    *m_logging_stream << "Step: " << step << " ";
    *m_logging_stream << "Movetype: " << movetype.m_label() << " ";
    *m_logging_stream << "Staples: " << m_origami_system.num_staples() << " ";
    *m_logging_stream << "Accepted: " << accepted << " ";
    *m_logging_stream << "Temp: " << m_origami_system.m_temp << " ";
    *m_logging_stream << "Energy: " << m_origami_system.energy() << " ";
//    *m_logging_stream << "Bias: " << m_system_bias.calc_bias() << " ";
    *m_logging_stream << "\n";
}

ConstantTGCMCSimulation::ConstantTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        GCMCSimulation(origami_system, params),
        m_steps {params.m_steps} {

    m_logging_stream = &cout;
    m_output_files = setup_output_files(params, params.m_output_filebase,
            m_origami_system);
}

AnnealingGCMCSimulation::AnnealingGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
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
        InputParameters& params) :
        GCMCSimulation(origami_system, params),
        m_num_reps {params.m_num_reps},
        m_exchange_interval {params.m_exchange_interval} {

    m_swaps = params.m_steps / m_exchange_interval;
    string string_rank {std::to_string(m_rank)};
    string output_filebase {params.m_output_filebase + "-" + string_rank};
    m_output_files = setup_output_files(params, output_filebase,
            m_origami_system);
    m_logging_stream = new ofstream {output_filebase + ".out"};

    // Initialize quantity index to replica index vectors
    if (m_rank == m_master_rep) {
        for (int i {0}; i != m_num_reps; i++) {
            m_q_to_repi.push_back(i);
        }
    }

    initialize_control_qs(params);

    // Initalize swap file
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
}

UTPTGCMCSimulation::UTPTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        PTGCMCSimulation(origami_system, params) {
    m_exchange_q_is.push_back(m_temp_i);
    m_exchange_q_is.push_back(m_staple_u_i);
}

HUTPTGCMCSimulation::HUTPTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        PTGCMCSimulation(origami_system, params) {
    m_exchange_q_is.push_back(m_temp_i);
    m_exchange_q_is.push_back(m_staple_u_i);
    m_exchange_q_is.push_back(m_bias_mult_i);
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

UmbrellaSamplingSimulation::UmbrellaSamplingSimulation(
        OrigamiSystem& origami,
        InputParameters& params) :
        GCMCSimulation(origami, params),
        m_params {params},
        m_num_iters {params.m_num_iters},
        m_steps {params.m_steps},
        // This is awfull
        m_system_order_params {dynamic_cast<OrigamiSystemWithBias*>(&origami)->
                get_system_order_params()},
        m_grid_bias {dynamic_cast<OrigamiSystemWithBias*>(&origami)->
            get_system_biases()->get_grid_bias()},
        m_max_D_bias {params.m_max_D_bias} {

    // Initalize log file
    m_solver_file.open(params.m_output_filebase + ".solver");

    // Setup grid bias
    //OrderParam* dist_sum {m_system_order_params->get_dist_sums()[0]};
    OrderParam* num_bound_domains {&m_system_order_params->get_num_bound_domains()};
    OrderParam* num_staples {&m_system_order_params->get_num_staples()};
    //m_grid_params = {dist_sum, num_bound_domains};
    //m_grid_params = {dist_sum, num_staples};
    m_grid_params = {num_bound_domains, num_staples};
    m_grid_bias->set_order_params(m_grid_params);

    // Set maximum allowed difference for determining if iteration is to be discarded
    m_equil_dif.push_back(5);
    m_equil_dif.push_back(2);
}

UmbrellaSamplingSimulation::~UmbrellaSamplingSimulation() {
}

void UmbrellaSamplingSimulation::run() {
    int n {0};

    // Production
    while (n != m_num_iters) {

        // Write each iteration's output to a seperate file
        string prefix {"-" + std::to_string(n)};
        string output_filebase {m_params.m_output_filebase + prefix};
        m_output_files = setup_output_files(m_params, output_filebase, m_origami_system);
        m_logging_stream = new ofstream {output_filebase + ".out"};

        simulate(m_steps);
        vector<GridPoint> new_points(m_s_i.size());
        vector<GridPoint> old_points(m_s_i.size());
        vector<GridPoint> old_only_points(m_S_n.size());
        fill_grid_sets(new_points, old_points, old_only_points);
        bool step_is_equil {iteration_equilibrium_step(new_points, old_points)};
        if (step_is_equil) {
            cout << "Discarding iteration\n\n";
            m_s_i.clear();
            m_f_i.clear();
            continue;
        }
        m_S_n.insert(m_s_i.begin(), m_s_i.end());
        update_grids(n, new_points, old_points, old_only_points);
        estimate_current_weights(n, new_points, old_only_points);
        //estimate_normalizations(n);
        update_bias();
        output_summary(n);
        n++;
        m_output_files.clear();
        delete m_logging_stream;
    }
}

void UmbrellaSamplingSimulation::update_grids(int n,
        vector<GridPoint> new_points, vector<GridPoint> old_points,
        vector<GridPoint> old_only_points) {
    // Update tracking stuff
    m_s_n.push_back(m_s_i);

    for (auto point: new_points) {
        vector<int> f_k_n(n);
        f_k_n.push_back(m_f_i[point]);
        m_f_n[point] = f_k_n;

        m_F_n[point] = f_k_n;

        vector<double> r_k_n(n);
        r_k_n.push_back(1);
        m_r_n[point] = r_k_n;

        vector<double> w_k_n(n);
        w_k_n.push_back(static_cast<double>(m_f_i[point]) / m_steps);
        m_w_n[point] = w_k_n;
    }
    for (auto point: old_points) {
        m_f_n[point].push_back(m_f_i[point]);

        int F_k_n {std::accumulate(m_f_n[point].begin(), m_f_n[point].end(), 0)};
        m_F_n[point].push_back(F_k_n);

        m_r_n[point].push_back(0);
        for (int i {0}; i != n + 1; i++) {
            m_r_n[point][i] = static_cast<double>(m_f_n[point][i]) / F_k_n;
        }

        m_w_n[point].push_back(static_cast<double>(m_f_i[point]) / m_steps);
    }

    for (auto point: old_only_points) {
        m_f_n[point].push_back(0);
        m_F_n[point].push_back(m_F_n[point].back());
        m_r_n[point].push_back(0);
        m_w_n[point].push_back(0);
    }

    m_s_i.clear();
    m_f_i.clear();
}

void UmbrellaSamplingSimulation::update_internal() {
    GridPoint point {};
    for (auto grid_param: m_grid_params) {
        point.push_back(grid_param->get_param());
    }
    m_s_i.insert(point);
    m_f_i[point] ++;
}

void UmbrellaSamplingSimulation::estimate_current_weights(int n,
        vector<GridPoint> new_points, vector<GridPoint> old_only_points) {
    // Average bias weights of iteration
    double ave_bias_weight {0};
    for (auto point: m_s_n.back()) {
        double point_bias {m_grid_bias->calc_bias(point)};
        ave_bias_weight += m_f_n[point].back() * std::exp(point_bias);
    }

    // Setup vectors for new points
    for (auto point: new_points) {
        vector<double> p_k_n(n);
        m_p_n[point] = p_k_n;
    }

    // Calculate for all visited points
    for (auto point: m_s_n.back()) {
        double point_bias {m_grid_bias->calc_bias(point)};
        double p_n_k {m_f_n[point].back() * std::exp(point_bias) /
            ave_bias_weight};
        m_p_n[point].push_back(p_n_k);
    }

    // Update vector for unvisted points
    for (auto point: old_only_points) {
        m_p_n[point].push_back(0);
    }
}

bool UmbrellaSamplingSimulation::iteration_equilibrium_step(
        vector<GridPoint> new_points, vector<GridPoint> old_points) {
    // Check if any visited gridpoint in current iteration further than x from
    // any previously visited gridpoints
    // Other mode advanced checks possible
    bool equilibrium_step {false};
    
    // This is a pretty ugly way of dealing with the first iteration step
    if (m_S_n.size() == 0) {
        equilibrium_step = false;
        return equilibrium_step;
    }

    // Must have at least one point that has been visited previously
    if (old_points.size() == 0) {
        equilibrium_step = true;
        return equilibrium_step;
    }

    // Check distance between new points and closest previously sampled points
    for (auto new_point: new_points) {
        // Does it matter what order I move through the dimensions?
        for (size_t dim {0}; dim != new_point.size(); dim++) {
            GridPoint closest_point {find_closest_point(m_S_n, new_point, dim)};
            if (abs(new_point[dim] - closest_point[dim]) > m_equil_dif[dim]) {
                equilibrium_step = true;
                break;
            }
        }
    }

    return equilibrium_step;
}

void UmbrellaSamplingSimulation::fill_grid_sets(vector<GridPoint>& new_points,
        vector<GridPoint>& old_points, vector<GridPoint>& old_only_points) {

    // First find states that were visited before and now
    vector<GridPoint>::iterator it;
    it = std::set_intersection(m_s_i.begin(), m_s_i.end(),
            m_S_n.begin(), m_S_n.end(), old_points.begin());
    old_points.resize(it - old_points.begin());

    // Now find newly sampled points
    it = std::set_difference(m_s_i.begin(), m_s_i.end(), old_points.begin(),
            old_points.end(), new_points.begin());
    new_points.resize(it - new_points.begin());

    // Now find points sampled before not in current
    it = std::set_difference(m_S_n.begin(), m_S_n.end(), old_points.begin(),
            old_points.end(), old_only_points.begin());
    old_only_points.resize(it - old_only_points.begin());

    return;
}

void UmbrellaSamplingSimulation::update_bias() {
    for (auto point: m_S_n) {

        // No T to be consistent with biases here
        double old_bias {m_E_w[point]};
        //double new_bias {std::log(m_P_n[point])};
        double p_k_n {m_p_n[point].back()};
        double D_bias;
        double new_bias;
        if (p_k_n == 0) {
            D_bias = -m_max_D_bias;
            new_bias = old_bias + D_bias;
        }
        else {
            new_bias = std::log(p_k_n);
            D_bias = new_bias - old_bias;
        }
        double updated_bias {new_bias};

        // Limit how quckly bias can change
        if (abs(D_bias) > m_max_D_bias) {
            if (D_bias > 0) {
                updated_bias = old_bias + m_max_D_bias;
            }
            else {
                updated_bias = old_bias - m_max_D_bias;
            }
        }

        m_E_w[point] = updated_bias;
    }
    m_grid_bias->replace_biases(m_E_w);
}

void UmbrellaSamplingSimulation::output_summary(int n) {
    cout << "Iteration: " << n << "\n";
    //cout << "Normalization constants:\n";
    //for (auto N: m_N) {
    //    cout << N << " ";
    //}
    cout << "\n";
    //cout << "Gridpoint w, r, p, P, E:\n";
    cout << "Gridpoint w, r, p, E:\n";
    for (auto point: m_S_n) {
        for (auto coor: point) {
            cout << coor << " ";
        }
        cout << std::setprecision(3);
        //cout << ": " << std::setw(10) << m_w_n[point][n] << std::setw(10) << m_r_n[point][n] << std::setw(10) << m_p_n[point][n] << std::setw(10) << m_P_n[point] << std::setw(10) << m_E_w[point] << "\n";
        cout << ": " << std::setw(10) << m_w_n[point][n] << std::setw(10) << m_r_n[point][n] << std::setw(10) << m_p_n[point][n] << std::setw(10) << m_E_w[point] << "\n";
    }
    cout << "\n";
}

void UmbrellaSamplingSimulation::estimate_normalizations(int n) {
    // Ugly special case
    if (n == 0) {
        m_N.push_back(1);
    }
    else {
        m_N.push_back(estimate_initial_normalization(n));
        estimate_final_normalizations();
    }

    // Update probabilities
    for (auto point: m_S_n) {
        double P_k_n {0};
        for (int i {0}; i != n + 1; i++) {
            P_k_n += m_r_n[point][i] * m_N[i] * m_p_n[point][i];
        }
        m_P_n[point] = P_k_n;
    }
}

double UmbrellaSamplingSimulation::estimate_initial_normalization(int n) {
    // Use a linear one-step optimization
    // Uses the previous Nis and estimates the current iteration's Ni
    
    // N_n = a(c)/b(c)

    // Calculate grid of c
    GridFloats c_grid {};
    long int F_prev_sum {n * m_steps};

    for (auto point: m_S_n) {
        double w_n_k {m_w_n[point][n]};
        double r_n_k {m_r_n[point][n]};
        double c {pow(r_n_k, 2) * m_F_n[point][n - 1] / F_prev_sum +
                w_n_k*pow((1 - r_n_k), 2)};
        c_grid[point] = c;
    }

    // Calculate a(c)
    double N_n_num {0};
    for (auto point: m_S_n) {
        N_n_num += m_P_n[point] * m_p_n[point][n] * c_grid[point];
    }

    // Calculate b(c)
    double N_n_denom {0};
    for (auto point: m_S_n) {
        N_n_denom += pow(m_p_n[point][n], 2) * c_grid[point];
    }

    return N_n_num / N_n_denom;
}

void UmbrellaSamplingSimulation::estimate_final_normalizations() {
    // With relative deviation square sum
    int num_iters {static_cast<int>(m_N.size())};
    double num_iters_as_double {static_cast<double>(m_N.size())};
    Problem problem;
    vector<int> parameter_block_sizes {num_iters - 1, num_iters, num_iters, 1, 1};
    for (int n {0}; n != num_iters; n++) {
        double n_as_double {static_cast<double>(n)};
        for (auto point: m_s_n[n]) {
            CostFunction* cost_function =
                    new RDevSquareSumCostFunction {1, parameter_block_sizes};

            /* Autodiff
            //DynamicAutoDiffCostFunction<RDevSquareSum, 4>* cost_function =
            //    new DynamicAutoDiffCostFunction<RDevSquareSum, 4>(
            //            new RDevSquareSum {});

            //cost_function->AddParameterBlock(num_iters);
            //cost_function->AddParameterBlock(num_iters);
            //cost_function->AddParameterBlock(num_iters);
            //cost_function->AddParameterBlock(1);
            //cost_function->AddParameterBlock(1);
            //cost_function->SetNumResiduals(1);
            */

            LossFunction* loss_function = new ScaledLoss {NULL, m_w_n[point][n],
                    ceres::TAKE_OWNERSHIP};
            // I am passing the core arrays that are held in vectors (&nvariable[0])
            problem.AddResidualBlock(cost_function, loss_function, &m_N[1],
                    &m_r_n[point][0], &m_p_n[point][0], &n_as_double, &num_iters_as_double);
            problem.SetParameterBlockConstant(&m_r_n[point][0]);
            problem.SetParameterBlockConstant(&m_p_n[point][0]);
            problem.SetParameterBlockConstant(&n_as_double);
            problem.SetParameterBlockConstant(&num_iters_as_double);
        }
        if (n != 0) {
            problem.SetParameterLowerBound(&m_N[1], n - 1, 0);
        }
    }

    Solver::Options options;
    options.max_num_iterations = 100;
    options.minimizer_type = ceres::TRUST_REGION;
    options.minimizer_progress_to_stdout = false;

    Solver::Summary summary;
    Solve(options, &problem, &summary);
    m_solver_file << "Iteration " << (num_iters - 1) << "\n";
    m_solver_file << summary.FullReport() << "\n\n";
}

RDevSquareSumCostFunction::RDevSquareSumCostFunction(int num_residuals,
        vector<int>& parameter_block_sizes) {
    set_num_residuals(num_residuals);
    vector<int>* parameter_block_sizes_ {mutable_parameter_block_sizes()};
    *parameter_block_sizes_ = parameter_block_sizes;
}

bool RDevSquareSumCostFunction::Evaluate(
        double const* const* params, // 0: N_n[1:]; 1: r_k_n; 2: p_k_n; 3: iter; 4 num_iters
        double* residuals,
        double** jacobians) const {

    // Calculate cost function
    double P_k_n {params[1][0] * params[2][0]};
    int num_iters {static_cast<int>(params[4][0])};
    for (int j {0}; j != num_iters - 1; j ++) {
        P_k_n += params[1][j + 1] * params[0][j] * params[2][j + 1];
    }
    
    int iter {static_cast<int>(params[3][0])};
    if (iter == 0) {
        residuals[0] = (params[2][0] - P_k_n) / P_k_n;
    }
    else {
        residuals[0] = (params[0][iter - 1]*params[2][iter] - P_k_n) / P_k_n;
    }

    // Calculate Jacobian of cost function
    if (jacobians != NULL && jacobians[0] != NULL) {
        double P_k_n_sqr {pow(P_k_n, 2)};
        if (iter == 0) {
            for (int j {0}; j != num_iters - 1; j++) {
                double term_two = -pow(params[2][iter], 2)*params[1][iter] /
                    P_k_n_sqr;
                jacobians[0][j] = term_two;
            }
        }
        else {
            for (int j {0}; j != num_iters - 1; j++) {
                double term_two = -params[0][iter - 1]*pow(params[2][iter], 2)*params[1][iter] /
                    P_k_n_sqr;
                if (j - 1 == iter) {
                    jacobians[0][j] = params[2][iter] / P_k_n + term_two;
                }
                else {
                    jacobians[0][j] = term_two;
                }
            }
        }
    }
    return true;
}

template <typename T>
bool RDevSquareSum::operator() (
        T const* const* params, // 0, N_n; 1, r_n_k; 2 p_n_k; 3, i; 4 num_iters
        T* residual) const {
    residual[0] = params[1][0] * params[0][0] * params[2][0];
    for (double j {1}; j != params[4][0]; j += 1) {
        int j_int {static_cast<int>(j)};
        residual[0] += params[1][j_int] * params[0][j_int] * params[2][j_int];
    }
    
    // Hack way of getting an index from the passed parameter
    double i;
    for (i = 0; i != params[3][0]; i++) {}
    int i_int = static_cast<int>(i);
    residual[0] = (params[0][i_int] * params[2][i_int] - residual[0]) / residual[0];

    return true;
}

GridPoint Simulation::find_closest_point(set<GridPoint> search_set,
        GridPoint target_point, int dim) {

    // I need a point from the search set to initialize things, but set has no
    // front method, so this is a hacky solution
    GridPoint closest_point;
    for (auto point: search_set) {
        closest_point = point;
        break;
    }
    int closest_dist {abs(target_point[dim] - closest_point[dim])};
    for (auto point: search_set) {
        int cur_dist {abs(target_point[dim] - point[dim])};
        if (cur_dist < closest_dist) {
            closest_point = point;
            closest_dist = cur_dist;
        }
    }

    return closest_point;
}
