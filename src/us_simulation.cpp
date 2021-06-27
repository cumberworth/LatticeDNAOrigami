// us_simulation.cpp

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
#include "boost/serialization/boost_unordered_map.hpp"
#include "boost/serialization/set.hpp"
#include "json/json.h"

#include "files.h"
#include "us_simulation.h"
#include "utility.h"

namespace us {

using std::ifstream;
using std::chrono::steady_clock;

using biasFunctions::SquareWellBiasFunction;
using files::OrigamiTrajInputFile;
using origami::Chains;

USGCMCSimulation::USGCMCSimulation(
        OrigamiSystem& origami,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        GCMCSimulation(origami, ops, biases, params),
        m_params {params},
        m_max_num_iters {params.m_max_num_iters},
        m_equil_steps {params.m_equil_steps},
        m_max_equil_dur {params.m_max_equil_dur},
        m_iter_steps {params.m_iter_steps},
        m_max_iter_dur {params.m_max_iter_dur},
        m_grid_bias {biases.get_grid_bias(params.m_us_grid_bias_tag)},
        m_max_D_bias {params.m_max_D_bias} {

    // Add reading in other internal variables here?
    // But I will also need to specify which iteration to begin on
    // Also I need to not clear the grid if I am loading it now

    // Read in weights if specified
    if (params.m_read_biases == true) {
        ifstream bias_file {params.m_biases_file};
        if (bias_file.good()) {
            *m_us_stream << "Reading biases from file\n";
            read_weights(params.m_biases_file);
            m_grid_bias.replace_biases(m_E_w);
        }
        else {
            *m_us_stream << "Bias file not good\n";
            throw utility::SimulationMisuse {};
        }
    }
    else {
        *m_us_stream << "No biases read in\n";
    }

    // Update starting configs if specified
    if (m_params.m_restart_from_config == true) {
        *m_us_stream << "Reading starting config from file\n";
        set_config_from_traj(
                m_params.m_restart_traj_file, params.m_restart_step);
    }
    else {
        *m_us_stream << "Starting from configuration in system file\n";
    }

    // Update internal grids if specified
    if (m_params.m_restart_us_iter == true) {
        *m_us_stream << "Restarting iteration from file\n";
        set_grids_from_file(m_params.m_restart_us_filebase);
    }
    else {
        *m_us_stream << "Starting new iteration\n";
    }
}

void USGCMCSimulation::run() {
    int n;
    run_equilibration();
    for (n = 0; n != m_max_num_iters; n++) {
        prepare_iteration(n);
        run_simulation(m_steps);
        process_iteration(n);
    }
}

void USGCMCSimulation::run_equilibration() {
    set_max_dur(m_max_equil_dur);

    // Setup output files
    string postfix {"_iter-equil"};
    string output_filebase {m_params.m_output_filebase + postfix};
    m_output_files = simulation::setup_output_files(
            m_params,
            output_filebase,
            m_origami_system,
            m_ops,
            m_biases,
            m_random_gens);
    m_logging_stream = new ofstream {output_filebase + ".out"};

    m_steps = m_equil_steps;
    m_steps = simulate(m_steps);

    // Cleanup
    clear_grids(); 
    close_output_files();
    delete m_logging_stream;
}

void USGCMCSimulation::prepare_iteration(int n) {
    set_max_dur(m_max_iter_dur);

    // Write each iteration's output to a seperate file
    string prefix {"_iter-" + std::to_string(n)};
    m_output_filebase = m_params.m_output_filebase + prefix;
    string filename {m_output_filebase + "-inp.biases"};
    output_weights(filename);
    m_output_files = simulation::setup_output_files(
            m_params,
            m_output_filebase,
            m_origami_system,
            m_ops,
            m_biases,
            m_random_gens);
    m_logging_stream = new ofstream {m_output_filebase + ".out"};

    // This will need to change for REMC
    m_steps = m_iter_steps;
}

void USGCMCSimulation::run_simulation(long long int steps) {

    // Need to change steps for REMC
    m_steps = simulate(steps);
}

void USGCMCSimulation::process_iteration(int n) {
    fill_grid_sets();
    m_S_n.insert(m_s_i.begin(), m_s_i.end());
    estimate_current_weights();
    update_grids(n);
    update_bias(n);
    output_summary(n);

    // Cleanup
    close_output_files();
    delete m_logging_stream;
    string filename {m_output_filebase + ".biases"};
    output_weights(filename);
    clear_grids();
}

void USGCMCSimulation::clear_grids() {
    m_s_i.clear();
    m_f_i.clear();
    m_w_i.clear();
    m_new_points.clear();
    m_old_points.clear();
    m_old_only_points.clear();
}

GridPoint USGCMCSimulation::get_current_point() {
    return m_grid_bias.get_point();
}

void USGCMCSimulation::read_weights(string filename) {
    ifstream jsonraw {filename, ifstream::binary};
    Json::Value jsonroot;
    jsonraw >> jsonroot;

    for (auto json_entry: jsonroot["biases"]) {
        GridPoint point {};
        for (auto comp: json_entry["point"]) {
            point.push_back(comp.asInt());
        }
        double point_bias {json_entry["bias"].asDouble()};
        m_E_w[point] = point_bias;
    }
}

void USGCMCSimulation::output_weights(string filename) {
    ofstream weights_file {filename};
    weights_file << "{\n";
    weights_file << "    \"biases\": [";
    for (auto iter = m_S_n.begin(); iter != m_S_n.end(); iter++) {

        // Last comma shit
        if (iter == m_S_n.begin()) {
            weights_file << "\n";
        }
        else {
            weights_file << ",\n";
        }

        GridPoint point {*iter};
        weights_file << "        {\n";
        weights_file << "            \"point\": [\n";
        for (size_t i {0}; i != point.size() - 1; i++) {
            weights_file << "                ";
            weights_file << std::to_string(point[i]);
            weights_file << ",\n";
        }
        weights_file << "                ";
        weights_file << std::to_string(point.back());
        weights_file << "\n";
        weights_file << "            ";
        weights_file << "],\n";
        weights_file << "            ";
        weights_file << "\"bias\": ";
        weights_file << std::to_string(m_E_w[point]);
        weights_file << "\n";
        weights_file << "        }";
    }
    weights_file << "\n    ]\n";
    weights_file << "}\n";
}

void USGCMCSimulation::set_config_from_traj(string filename, int step) {
    OrigamiTrajInputFile traj_inp {filename};
    Chains config {traj_inp.read_config(step)};
    m_origami_system.set_config(config);
}

void USGCMCSimulation::set_config_from_chains(Chains chains) {
    m_origami_system.set_config(chains);
}

void USGCMCSimulation::set_output_stream(ostream* out_stream) {
    m_us_stream = out_stream;
}

int USGCMCSimulation::get_grid_dim() { return m_grid_bias.get_dim(); }

void USGCMCSimulation::update_internal(long long int step) {
    GridPoint point {m_grid_bias.get_point()};
    m_s_i.insert(point);
    m_f_i[point]++;

    // Write archive of US state
    // Probably shouldn't use this output freq
    if (m_params.m_order_params_output_freq != 0 and
        step % m_params.m_order_params_output_freq == 0) {

        std::ofstream S_n_outfile {m_output_filebase + ".S_n"};
        boost::archive::text_oarchive S_n_archive {S_n_outfile};
        S_n_archive << m_S_n;

        std::ofstream s_i_outfile {m_output_filebase + ".s_i"};
        boost::archive::text_oarchive s_i_archive {s_i_outfile};
        s_i_archive << m_s_i;

        std::ofstream f_i_outfile {m_output_filebase + ".f_i"};
        boost::archive::text_oarchive f_i_archive {f_i_outfile};
        f_i_archive << m_f_i;
    }
}

void USGCMCSimulation::estimate_current_weights() {
    // Average bias weights of iteration
    double ave_bias_weight {0};
    for (auto point: m_s_i) {
        double point_bias {m_grid_bias.calc_bias(point)};
        ave_bias_weight += m_f_i[point] * std::exp(point_bias);
    }

    // Calculate for all visited points
    for (auto point: m_s_i) {
        double point_bias {m_grid_bias.calc_bias(point)};
        double p_n_k {m_f_i[point] * std::exp(point_bias) / ave_bias_weight};
        m_p_i[point] = p_n_k;
    }

    // Update vector for unvisted points
    for (auto point: m_old_only_points) {
        m_p_i[point] = 0;
    }
}

void USGCMCSimulation::set_grids_from_file(string filebase) {
    std::ifstream S_n_infile {filebase + ".S_n"};
    boost::archive::text_iarchive S_n_archive {S_n_infile};
    S_n_archive >> m_S_n;

    std::ifstream s_i_infile {filebase + ".s_i"};
    boost::archive::text_iarchive s_i_archive {s_i_infile};
    s_i_archive >> m_s_i;

    std::ifstream f_i_infile {filebase + ".f_i"};
    boost::archive::text_iarchive f_i_archive {f_i_infile};
    f_i_archive >> m_f_i;
    for (auto f: m_f_i) {
        *m_us_stream << f.first[0] << ": " << f.second << std::endl;
    }
}

void USGCMCSimulation::fill_grid_sets() {
    m_new_points.resize(m_s_i.size());
    m_old_points.resize(m_s_i.size());
    m_old_only_points.resize(m_S_n.size());

    // First find states that were visited before and now
    vector<GridPoint>::iterator it;
    it = std::set_intersection(
            m_s_i.begin(),
            m_s_i.end(),
            m_S_n.begin(),
            m_S_n.end(),
            m_old_points.begin());
    m_old_points.resize(it - m_old_points.begin());

    // Now find newly sampled points
    it = std::set_difference(
            m_s_i.begin(),
            m_s_i.end(),
            m_old_points.begin(),
            m_old_points.end(),
            m_new_points.begin());
    m_new_points.resize(it - m_new_points.begin());

    // Now find points sampled before not in current
    it = std::set_difference(
            m_S_n.begin(),
            m_S_n.end(),
            m_old_points.begin(),
            m_old_points.end(),
            m_old_only_points.begin());
    m_old_only_points.resize(it - m_old_only_points.begin());

    return;
}

SimpleUSGCMCSimulation::SimpleUSGCMCSimulation(
        OrigamiSystem& origami,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        USGCMCSimulation(origami, ops, biases, params) {}

void SimpleUSGCMCSimulation::update_bias(int) {
    for (auto point: m_S_n) {

        // No T to be consistent with biases here
        double old_bias {m_E_w[point]};
        double p_k_n {m_p_i[point]};
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
    m_grid_bias.replace_biases(m_E_w);
}

void SimpleUSGCMCSimulation::update_grids(int) {
    for (auto point: m_s_i) {
        m_w_i[point] = (static_cast<double>(m_f_i[point]) / m_steps);
    }
    for (auto point: m_old_only_points) {
        m_w_i[point] = 0;
    }
}

void SimpleUSGCMCSimulation::output_summary(int n) {
    *m_us_stream << "Iteration: " << n << "\n";
    *m_us_stream << "\n";
    *m_us_stream << "Gridpoint w, P, E:\n";
    for (auto point: m_S_n) {
        for (auto coor: point) {
            *m_us_stream << coor << " ";
        }
        *m_us_stream << std::setprecision(3);
        *m_us_stream << ": ";
        *m_us_stream << std::setw(10) << m_w_i[point];
        *m_us_stream << std::setw(10) << m_p_i[point];
        *m_us_stream << std::setw(10) << m_E_w[point];
        *m_us_stream << "\n";
    }
    *m_us_stream << "\n";
}

MWUSGCMCSimulation::MWUSGCMCSimulation(
        OrigamiSystem& origami,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        GCMCSimulation {origami, ops, biases, params},
        m_params {params},
        m_max_num_iters {params.m_max_num_iters},
        m_output_filebase {params.m_output_filebase} {

    parse_windows_file(params.m_windows_file);
    setup_window_variables();
    setup_window_restraints();
    setup_window_sims(origami);
}

MWUSGCMCSimulation::~MWUSGCMCSimulation() {
    delete m_us_sim;
    delete m_us_stream;
}

void MWUSGCMCSimulation::run() {
    int n;
    m_us_sim->run_equilibration();
    for (n = 0; n != m_max_num_iters; n++) {
        m_us_sim->prepare_iteration(n);
        m_us_sim->run_simulation(m_us_sim->m_steps);
        m_us_sim->process_iteration(n);
        if (m_rank == m_master_node) {
            output_iter_summary(n);
        }
    }
}

void MWUSGCMCSimulation::setup_window_variables() {
    for (int i {0}; i != m_windows; i++) {
        string window_postfix {"_win-"};
        for (auto j: m_window_mins[i]) {
            window_postfix += std::to_string(j);
            window_postfix += "-";
        }
        for (auto j: m_window_maxs[i]) {
            window_postfix += "-";
            window_postfix += std::to_string(j);
        }
        m_window_postfixes.push_back(window_postfix);
        string output_filebase {m_params.m_output_filebase + window_postfix};
        m_output_filebases.push_back(output_filebase);
    }
}

void MWUSGCMCSimulation::setup_window_restraints() {
    GridPoint window_min {m_window_mins[m_rank]};
    GridPoint window_max {m_window_maxs[m_rank]};

    for (size_t i {0}; i != m_window_bias_tags.size(); i++) {
        string bias_tag {m_window_bias_tags[i]};
        SquareWellBiasFunction& bias_f {
                m_biases.get_square_well_bias(bias_tag)};
        int op_min {window_min[i]};
        int op_max {window_max[i]};
        bias_f.set_min_op(op_min);
        bias_f.set_max_op(op_max);
    }
}

void MWUSGCMCSimulation::setup_window_sims(OrigamiSystem& origami) {

    // US sims reads filebase names from params object, so update
    string window_postfix {m_window_postfixes[m_rank]};
    string output_filebase {m_output_filebases[m_rank]};
    m_params.m_output_filebase = output_filebase;

    // If available read modify filename for input biases
    if (m_params.m_read_biases == true) {
        m_params.m_biases_file = m_params.m_biases_filebase + window_postfix;
        m_params.m_biases_file += ".biases";
    }

    // If available modify filename for seperate starting config for each window
    // Standard names
    if (m_params.m_restart_from_config == true) {
        m_params.m_restart_traj_filebase += window_postfix;
        m_params.m_restart_traj_file = m_params.m_restart_traj_filebase +
                                       m_params.m_restart_traj_postfix;
    }
    if (m_params.m_restart_us_iter == true) {
        m_params.m_restart_us_filebase += window_postfix;
    }

    // Names invididually specified in param file
    if (not m_params.m_restart_traj_files.empty()) {
        m_params.m_restart_traj_file = m_params.m_restart_traj_files[m_rank];
        m_params.m_restart_step = m_params.m_restart_steps[m_rank];
    }

    // Create simulation objects
    m_us_sim = new SimpleUSGCMCSimulation {origami, m_ops, m_biases, m_params};
    m_us_stream = new ofstream {output_filebase + ".out"};
    m_us_sim->set_output_stream(m_us_stream);
    m_grid_dim = m_us_sim->get_grid_dim();
}

void MWUSGCMCSimulation::output_iter_summary(int n) {
    cout << "Iteration: " << n << "\n";
}

void MWUSGCMCSimulation::parse_windows_file(string filename) {
    ifstream file {filename};
    string window_raw;

    // Get op tags
    string window_bias_tags_raw;
    std::getline(file, window_bias_tags_raw);
    std::istringstream tags_stream {window_bias_tags_raw};
    while (not tags_stream.eof()) {
        string tag;
        tags_stream >> tag;
        m_window_bias_tags.push_back(tag);
    }

    // Get mins and maxes
    while (std::getline(file, window_raw)) {
        std::istringstream window_stream {window_raw};
        string window_min_raw;
        std::getline(window_stream, window_min_raw, ',');
        std::istringstream window_min_stream {window_min_raw};
        GridPoint min_point {};
        while (not window_min_stream.eof()) {
            int comp;
            window_min_stream >> comp;
            min_point.push_back(comp);
        }
        m_window_mins.push_back(min_point);

        GridPoint max_point {};
        while (not window_stream.eof()) {
            int comp;
            window_stream >> comp;
            max_point.push_back(comp);
        }
        m_window_maxs.push_back(max_point);

        m_windows++;
    }
}

PTMWUSGCMCSimulation::PTMWUSGCMCSimulation(
        OrigamiSystem& origami,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        MWUSGCMCSimulation {origami, ops, biases, params},
        m_iter_swaps {params.m_iter_swaps},
        m_max_iter_dur {params.m_max_iter_dur},
        m_exchange_interval {params.m_exchange_interval},
        m_attempt_count(m_windows - 1, 0),
        m_swap_count(m_windows - 1, 0),
        m_win_to_configi(m_windows, 0) {
    initialize_swap_file(m_params);
}

void PTMWUSGCMCSimulation::run() {

    // Somehow based on input file decide where to restart
    int n;
    if (m_params.m_restart_us_iter == false) {
        m_us_sim->run_equilibration();
    }
    for (n = 0; n != m_max_num_iters; n++) {
        m_us_sim->prepare_iteration(n);
        run_swaps(m_iter_swaps, m_max_iter_dur);
        m_us_sim->process_iteration(n);
        if (m_rank == m_master_node) {
            output_iter_summary(n);
        }
        std::fill(m_attempt_count.begin(), m_attempt_count.end(), 0);
        std::fill(m_swap_count.begin(), m_swap_count.end(), 0);
    }
}

void PTMWUSGCMCSimulation::run_swaps(long long int swaps, long long int dur) {
    m_us_sim->m_steps = 0;
    auto start = steady_clock::now();
    for (int swap_i {1}; swap_i != m_iter_swaps + 1; swap_i++) {
        m_us_sim->simulate(m_exchange_interval, m_us_sim->m_steps, false);
        m_us_sim->m_steps += m_exchange_interval;

        // This is different than other time kills as will keep steps past the
        // duration as things are already save internally
        if (m_rank == m_master_node) {
            std::chrono::duration<double> dt {(steady_clock::now() - start)};
            if (dt.count() > dur) {
                master_send_kill(swap_i);
                cout << "Maximum time allowed reached\n";
                break;
            }
        }

        if (m_rank != m_master_node) {
            slave_send_ops(swap_i);
            if (not slave_send_and_recieve_chains(swap_i)) {
                break;
            }
        }
        else {
            write_swap_entry(m_us_sim->m_steps);
            attempt_exchange(swap_i);
        }
    }

    // Write end-of-simulation data
    if (m_rank == m_master_node) {
        write_swap_entry(m_us_sim->m_steps);
        write_acceptance_freqs();
        m_swapfile.close();
    }
}

void PTMWUSGCMCSimulation::slave_send_ops(int swap_i) {
    GridPoint point {m_us_sim->get_current_point()};
    cout << "Win " << m_rank << ": Sending point to " << m_master_node << "\n";
    m_world.send(m_master_node, swap_i, point);
    cout << "Win " << m_rank << ": Sent point to " << m_master_node << "\n";
    cout << "Win " << m_rank << ": Sending biases to " << m_master_node << "\n";
    m_world.send(m_master_node, swap_i, m_us_sim->m_E_w);
    cout << "Win " << m_rank << ": Sent biases to " << m_master_node << "\n";
}

bool PTMWUSGCMCSimulation::slave_send_and_recieve_chains(int swap_i) {
    int win_i;
    cout << "Win " << m_rank << ": Recieving win_i from " << m_master_node
         << "\n";
    m_world.recv(m_master_node, swap_i, win_i);
    cout << "Win " << m_rank << ": Recieved win_i from " << m_master_node
         << "\n";
    if (win_i == 999) {
        return false;
    }
    else if (win_i > m_rank) {
        Chains chains_send {m_us_sim->get_chains()};
        cout << "Win " << m_rank << ": Sending chains (size "
             << chains_send.size() << ") to " << win_i << "\n";
        m_world.send(win_i, swap_i, chains_send);
        cout << "Win " << m_rank << ": Sent chains to " << win_i << "\n";
        cout << "Win " << m_rank << ": Recieving chains from " << win_i << "\n";
        Chains chains_rec;
        m_world.recv(win_i, swap_i, chains_rec);
        cout << "Win " << m_rank << ": Recieved chains (size "
             << chains_rec.size() << ") from " << win_i << "\n";
        m_us_sim->set_config_from_chains(chains_rec);
    }
    else if (win_i < m_rank) {
        Chains chains_rec;
        cout << "Win " << m_rank << ": Recieving chains from " << win_i << "\n";
        m_world.recv(win_i, swap_i, chains_rec);
        cout << "Win " << m_rank << ": Recieved chains (size "
             << chains_rec.size() << ") from " << win_i << "\n";
        Chains chains_send {m_us_sim->get_chains()};
        cout << "Win " << m_rank << ": Sending chains (size "
             << chains_send.size() << ") to " << win_i << "\n";
        m_world.send(win_i, swap_i, chains_send);
        cout << "Win " << m_rank << ": Sent chains to " << win_i << "\n";
        m_us_sim->set_config_from_chains(chains_rec);
    }

    return true;
}

void PTMWUSGCMCSimulation::master_send_kill(int swap_i) {
    for (int i {1}; i != m_windows; i++) {
        cout << "Win " << m_rank << ": About to send kill signal to " << i
             << "\n";
        m_world.send(i, swap_i, 999);
    }
}

void PTMWUSGCMCSimulation::attempt_exchange(int swap_i) {
    cout << "Swap " << swap_i << "\n";

    // Collect order params
    vector<GridPoint> points {m_us_sim->get_current_point()};
    vector<GridFloats> biases {m_us_sim->m_E_w};
    vector<int> win_to_win {0};
    for (int i {1}; i != m_windows; i++) {
        win_to_win.push_back(i);
        GridPoint point;
        cout << "Win " << m_rank << ": Recieving point from " << i << "\n";
        m_world.recv(i, swap_i, point);
        cout << "Win " << m_rank << ": Recieved point from " << i << "\n";
        points.push_back(point);
        GridFloats bias;
        cout << "Win " << m_rank << ": Recieving biases from " << i << "\n";
        m_world.recv(i, swap_i, bias);
        cout << "Win " << m_rank << ": Recieved biases from " << i << "\n";
        biases.push_back(bias);
    }

    // Iterate through pairs in current set and attempt swap
    int swap_set {swap_i % 2};
    for (int i {swap_set}; i < (m_windows - 1); i += 2) {
        m_attempt_count[i]++;
        int win_1 {i};
        int win_2 {i + 1};
        GridPoint point_1 {points[win_1]};
        GridPoint point_2 {points[win_2]};

        // Test if both inside others window
        bool inside_windows {true};
        for (size_t j {0}; j != m_grid_dim; j++) {
            int min_comp {m_window_mins[win_1][j]};
            int max_comp {m_window_maxs[win_1][j]};
            if (point_2[j] < min_comp or point_2[j] > max_comp) {
                inside_windows = false;
                break;
            }
            min_comp = m_window_mins[win_2][j];
            max_comp = m_window_maxs[win_2][j];
            if (point_1[j] < min_comp or point_1[j] > max_comp) {
                inside_windows = false;
                break;
            }
        }

        if (inside_windows) {
            bool accepted {true};
            if (point_1 != point_2) {
                GridFloats biases_1 {biases[win_1]};
                GridFloats biases_2 {biases[win_2]};
                double point_1_bias_1 {biases_1[point_1]};
                double point_2_bias_1 {biases_1[point_2]};
                double point_1_bias_2 {biases_2[point_1]};
                double point_2_bias_2 {biases_2[point_2]};
                double point_1_bias_diff {point_1_bias_1 - point_1_bias_2};
                double point_2_bias_diff {point_2_bias_2 - point_2_bias_1};
                double point_bias_diff_sum {
                        point_1_bias_diff + point_2_bias_diff};
                double p_accept {std::min({1.0, exp(point_bias_diff_sum)})};
                accepted = test_acceptance(p_accept);
            }
            if (accepted) {
                m_swap_count[i]++;
                win_to_win[win_1] = win_2;
                win_to_win[win_2] = win_1;
            }
        }
    }
    for (int i {1}; i != m_windows; i++) {
        cout << "Win " << m_rank << ": Sending win_i to " << i << "\n";
        m_world.send(i, swap_i, win_to_win[i]);
        cout << "Win " << m_rank << ": Sent win_i to " << i << "\n";
    }
    if (win_to_win[0] != 0) {
        Chains chains_send {m_us_sim->get_chains()};
        cout << "Win " << m_rank << ": Sending chains (size "
             << chains_send.size() << ") to " << win_to_win[0] << "\n";
        m_world.send(win_to_win[0], swap_i, chains_send);
        cout << "Win " << m_rank << ": Sent chains to " << win_to_win[0]
             << "\n";
        Chains chains_rec;
        cout << "Win " << m_rank << ": Recieving chains from " << win_to_win[0]
             << "\n";
        m_world.recv(win_to_win[0], swap_i, chains_rec);
        cout << "Win " << m_rank << ": Recieved chains (size "
             << chains_rec.size() << ") from " << win_to_win[0] << "\n";
        m_us_sim->set_config_from_chains(chains_rec);
    }
    m_win_to_configi = win_to_win;
}

bool PTMWUSGCMCSimulation::test_acceptance(double p_accept) {
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

void PTMWUSGCMCSimulation::write_acceptance_freqs() {

    for (size_t i {0}; i != m_attempt_count.size(); i++) {
        cout << i << " ";
        cout << i + 1 << " ";
        cout << m_swap_count[i] << " ";
        cout << m_attempt_count[i] << " ";
        cout << (static_cast<double>(m_swap_count[i]) / m_attempt_count[i])
             << " ";
        cout << "\n";
    }
    cout << "\n";
}

void PTMWUSGCMCSimulation::initialize_swap_file(InputParameters& params) {
    if (m_rank == m_master_node) {
        m_swapfile.open(m_output_filebase + ".swp");
        for (int rep {0}; rep != m_windows; rep++) {
            m_swapfile << rep << " ";
            m_win_to_configi[rep] = rep;
        }
        m_swapfile << "\n";
    }
}

void PTMWUSGCMCSimulation::write_swap_entry(long long int step) {
    if (step % m_params.m_configs_output_freq == 0) {
        for (auto repi: m_win_to_configi) {
            m_swapfile << repi << " ";
        }
        m_swapfile << "\n";
    }
}
} // namespace us
