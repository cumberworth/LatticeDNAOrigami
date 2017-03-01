// simulation.cpp

#include <random>
#include <iostream>
#include <set>
#include <algorithm>

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
        m_bias_mults = params.m_bias_mults;
        for (int i {0}; i != m_num_reps; i++) {
            m_q_to_repi.push_back(i);
        }
    }

    // Initalize swap file
    if (m_rank == m_master_rep) {
        m_swapfile.open(params.m_output_filebase + ".swp");
        for (int rep {0}; rep != m_num_reps; rep++) {
            for (size_t q {0}; q != m_control_qs.size(); q++) {
                m_swapfile <<  m_control_qs[q][rep];
                if (q != m_control_qs.size()) {
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

UTPTGCMCSimulation::UTPTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        PTGCMCSimulation(origami_system, params) {
    if (m_rank == m_master_rep) {
        m_control_qs.push_back(params.m_temps);
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
    }

    // Initialize quantities of each replica (updating on origami happens in run)
    for (int i {0}; i != m_num_reps; i++) {
        if (m_rank == i) {
            m_temp = params.m_temps[i];
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
}

//            m_bias_mult = params.m_bias_mults[i];

PTGCMCSimulation::~PTGCMCSimulation() {
    delete m_logging_stream;
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
            slave_recieve(swap_i);
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

void TPTGCMCSimulation::update_control_qs() {
        m_origami_system.update_temp(m_temp);
}

void UTPTGCMCSimulation::update_control_qs() {
        m_origami_system.update_temp(m_temp);
        m_origami_system.update_staple_u(m_staple_u);
}
        //m_system_bias.update_bias_mult(m_bias_mult);

void TPTGCMCSimulation::slave_send(int swap_i) {
    // Send quantities to master
    ThermoOfHybrid DH_DS {m_origami_system.enthalpy_and_entropy()};
    m_world.send(m_master_rep, swap_i, DH_DS.enthalpy);
}


void TPTGCMCSimulation::slave_recieve(int swap_i) {
    // Receive quantities to master
    m_world.recv(m_master_rep, swap_i, m_temp);
}

void UTPTGCMCSimulation::slave_send(int swap_i) {
    // Send quantities to master
    ThermoOfHybrid DH_DS {m_origami_system.enthalpy_and_entropy()};
    int N {m_origami_system.num_staples()};

    // Note the order this is sent is important, consider making this part of
    // the code by iterating through something that both slave send and master
    // receieve have
    m_world.send(m_master_rep, swap_i, DH_DS.enthalpy);
    m_world.send(m_master_rep, swap_i, N);
}

void UTPTGCMCSimulation::slave_recieve(int swap_i) {
    // Receive quantities to master
    m_world.recv(m_master_rep, swap_i, m_temp);
    m_world.recv(m_master_rep, swap_i, m_staple_u);
}

    //double bias {m_system_bias.calc_bias()};
    //m_world.send(m_master_rep, swap_i, bias);
    //m_world.recv(m_master_rep, swap_i, m_bias_mult);

void PTGCMCSimulation::master_receive(int swap_i,
        vector<vector<double>>& dependent_qs) {
    master_get_dependent_qs(dependent_qs);
    // Receive dependent quantities from slaves
    for (int rep_i {1}; rep_i != m_num_reps; rep_i++) {
        for (auto dependent_q: dependent_qs) {
            double q;
            m_world.recv(rep_i, swap_i, q);
            dependent_q.push_back(q);
        }
    }
}

void PTGCMCSimulation::master_send(int swap_i) {
    // Send temperatures and chemical potentials to slaves
    for (int q_i {0}; q_i != m_num_reps; q_i++) {
        int rep_i {m_q_to_repi[q_i]};
        if (rep_i == m_master_rep) {
                master_set_control_qs();
        }
        else {
            for (auto control_q: m_control_qs) {
                m_world.send(rep_i, swap_i, control_q[q_i]);
            }
        }
    }
}

void TPTGCMCSimulation::master_set_control_qs() {
    // I could eliminate this by having the replica variables also replaced with a single vector
    m_temp = m_control_qs[m_temp_i][m_master_rep];
}

void UTPTGCMCSimulation::master_set_control_qs() {
    // I could eliminate this by having the replica variables also replaced with a single vector
    m_temp = m_control_qs[m_temp_i][m_master_rep];
    m_staple_u = m_control_qs[m_staple_u_i][m_master_rep];
}

void PTGCMCSimulation::master_get_dependent_qs(
        vector<vector<double>>& dependent_qs) {
    ThermoOfHybrid DH_DS {m_origami_system.enthalpy_and_entropy()};
    double staples {static_cast<double>(m_origami_system.num_staples())};
    dependent_qs[0].push_back(DH_DS.enthalpy);
    dependent_qs[1].push_back(staples);
}

    //vector<double> biases {m_system_bias.calc_bias()};

void PTGCMCSimulation::attempt_exchange(int swap_i,
        vector<int>& attempt_count, vector<int>& swap_count) {
    // Alternates between two sets of pairs, allows for changes in T, u, and N

    // Collect results from all replicas
    vector<vector<double>> dependent_qs {};
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
        for (auto control_q: m_control_qs) {
            double q_1 {control_q[repi1]};
            double q_2 {control_q[repi2]};
            dependent_q_pairs.push_back({q_1, q_2});
        }

        double p_accept {calc_acceptance_p(control_q_pairs, dependent_q_pairs)};
        bool accept {test_acceptance(p_accept)};

        // If accepted swap indices from quantities to replicas
        if (accept) {
            swap_count[i]++;
            int repi1_old {repi1};
            repi1 = repi2;
            repi2 = repi1_old;
            m_q_to_repi[i] = repi1;
            m_q_to_repi[i + 1] = repi2;
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

double UTPTGCMCSimulation::calc_acceptance_p(
        vector<pair<double, double>> control_q_pairs,
        vector<pair<double, double>> dependent_q_pairs) {

    double temp1 {control_q_pairs[m_temp_i].first};
    double temp2 {control_q_pairs[m_temp_i].second};
    double staple_u1 {control_q_pairs[m_temp_i].first};
    double staple_u2 {control_q_pairs[m_temp_i].second};
    double enthalpy1 {dependent_q_pairs[m_enthalpy_i].first};
    double enthalpy2 {dependent_q_pairs[m_enthalpy_i].second};
    double N1 {dependent_q_pairs[m_staples_i].first};
    double N2 {dependent_q_pairs[m_staples_i].second};

    // Energies are actually E/B, so multiply by T
    double DB {1/temp2 - 1/temp1};
    double DH {enthalpy2 - enthalpy1};
    double DN {N2 - N1};
    double DBU {staple_u2 / temp2 - staple_u1 / temp1};
    double p_accept {min({1.0, exp(DB*(DH) - DBU*DN)})};

    return p_accept;
}

    //double DBias {bias2 - bias1};

void PTGCMCSimulation::write_swap_entry() {
    for (auto repi: m_q_to_repi) {
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

UmbrellaSamplingSimulation::UmbrellaSamplingSimulation(OrigamiSystem& origami,
        InputParameters& params) :
        GCMCSimulation(origami, params),
        m_system_order_params {params, origami},
        m_system_biases {origami, m_system_order_params, params} {

    // Initalize ceres solver file
    m_solver_file.open(params.m_output_filebase + ".solver");
}

UmbrellaSamplingSimulation::~UmbrellaSamplingSimulation() {
}

void UmbrellaSamplingSimulation::run() {
    int n {0};
    while (n != m_num_iters) {
        simulate(m_steps_per_iter);
        estimate_current_weights();
        if (iteration_equilibrium_step()) {
            // what do with internal variables?
            continue;
        }
        estimate_normalizations(n);
        update_biases();
        n++;
    }
}

void UmbrellaSamplingSimulation::update_internal() {
    GridPoint point {something with the order params};
    // do i want to put it into an isolated first?
    m_s_n.back().insert(point);
    m_f_n[point].back() ++;
}

void UmbrellaSamplingSimulation::estimate_current_weights() {
    // Average bias weights of iteration
    double sum_of_biases {0};
    for (auto point: m_s_n.back()) {
        sum_of_biases += ?;
    }
    for (auto point: m_s_n.back()) {
        double p_n_k {m_f_n[point].back() * ? / sum_of_biases};
        m_p_n[point].push_back(p_n_k);
    }
}

bool UmbrellaSamplingSimulation::iteration_equilibrium_step() {
    // Check if any visited gridpoint in current iteration further than x from
    // any previously visited gridpoints
    // Other mode advanced checks possible
    bool equilibrium_step {false};

    // First find intersection of previously and currently visited states
    vector<GridPoint> it {};
    std::set_intersection(m_s_n.back().begin(), m_s_n.back().end(),
            m_S_n.begin(), m_S_n.end(), it.front());

    // Now find newly sampled points
    vector<GridPoint> dif {};
    std::set_difference(it.begin(), it.end(), m_s_n.back().begin(),
            m_s_n.back().end(), dif.front());

    for (int d {0}; d != dif.size(); d++) {
        int closet_point {find_closest_point(m_S_n, d)};
        if (abs(dif[d] - closest_point[d]) > m_equil_dif[d]) {
            equilibrium_step = true;
            break;
        }

    return equilibrium_step;
}

void UmbrellaSamplingSimulation::estimate_normalizations(int n) {
    estimate_initial_normalization(int n);
    estimate_final_normalizations();
}

double UmbrellaSamplingSimulation::estimate_initial_normalization(int n) {
    // Use a linear one-step optimization
    // Uses the previous Nis and estimates the current iteration's Ni
    
    // N_n = a(c)/b(c)

    // Calculate grid of c
    GridFloats c_grid {};
    int F_prev_sum {0};
    for (auto point: m_S_n) {
        F_prev_sum += m_F_n[point];
    }

    for (auto point: m_S_n) {
        double w_n_k {m_w_n[point][n]};
        double r_n_k {m_r_n[point][n]};
        c_grid[point] = m_F_n[point] / pow(F_prev_sum*r_n_k, 2) +
                w_n_k*pow((1 - r_n_k), 2);
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
    for (int n {0}; n != num_iters; n++) {
        for (auto point: m_s_n[n]) {
            DynamicAutoDiffCostFunction<RDevSquareSum, 4>* cost_function =
                    new DynamicAutoDiffCostFunction<RDevSquareSum, 4> {
                            new RDevSquareSum {}};

            cost_function->AddParameterBlock(num_iters);
            cost_function->AddParameterBlock(num_iters);
            cost_function->AddParameterBlock(num_iters);
            cost_function->AddParameterBlock(1);

            LossFunction* loss_function = new ScaledLoss {NULL, m_w_n[point][n],
                    ceres::TAKE_OWNERSHIP};
            problem.AddResidualBlock(cost_function, loss_function, &m_N[0],
                    &m_r_n[point][0], &m_p_n[point][0], &m_N[n], &m_p_n[point][n],
                    &num_iters_as_double);
            problem.SetParameterBlockConstant(&m_N[0]);
            problem.SetParameterBlockConstant(&m_r_n[point][0]);
            problem.SetParameterBlockConstant(&m_p_n[point][0]);
            problem.SetParameterBlockConstant(&m_N[n]);
            problem.SetParameterBlockConstant(&m_p_n[point][n]);
            problem.SetParameterBlockConstant(&num_iters_as_double);
        }
    }

    Solver::Options options;
    options.max_num_iterations = 100;
    options.minimizer_type = ceres::LINE_SEARCH;
    options.minimizer_progress_to_stdout = false;

    Solver::Summary summary;
    Solve(options, &problem, &summary);
    m_solver_file << "Iteration " << (num_iters - 1) << "\n";
    m_solver_file << summary.FullReport() << "\n\n";
}

template <typename T>
bool RDevSquareSum::operator() (
        T const* const* parameters, // 0, N_n; 1, r_n_k; 2 p_n_k; 3, N_i; 4, p_i_k; 5, num_iters
        T* residual) const {
    residual[0] = parameters[1][0] * parameters[0][0] * parameters[2][0];
    for (double j {1}; j != parameters[5][0]; j += 1) {
        int j_int {static_cast<int>(j)};
        residual[0] += parameters[1][j_int] * parameters[0][j_int] * parameters[2][j_int];
    }
    residual[0] = parameters[4][0] * parameters[5][0] - residual[0];

    return true;
    GridOfInts a {};
}

Simulation::find_closest_point(set<GridPoint> search_set,
        GridPoint target_point, int dim) {
    GridPoint closest_point {search_set.front()};
    int closest_dist {abs(target_point[d] - closest_point[d])};
    for (auto point: search_set) {
        int cur_dist {abs(target_point[d] - point[d])};
        if (cur_dist < closest_dist) {
            closest_point = point;
            closest_dist = cur_dist;
        }
    }

    return closest_point;
}
