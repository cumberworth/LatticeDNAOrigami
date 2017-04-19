// us_simulation.cpp

#include <cmath>

#include "json/json.h"
#include "us_simulation.h"

using namespace US;
using namespace Simulation;

GridPoint US::find_closest_point(set<GridPoint> search_set,
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

USGCMCSimulation::USGCMCSimulation(
        OrigamiSystem& origami,
        InputParameters& params) :
        GCMCSimulation(origami, params),
        m_params {params},
        m_max_num_iters {params.m_max_num_iters},
        m_equil_steps {params.m_equil_steps},
        m_steps {params.m_steps},
        m_prod_steps {params.m_prod_steps},
        m_max_rel_P_diff {params.m_max_rel_P_diff},
        // This is awfull
        m_system_order_params {dynamic_cast<OrigamiSystemWithBias*>(&origami)->
                get_system_order_params()},
        m_grid_bias {dynamic_cast<OrigamiSystemWithBias*>(&origami)->
            get_system_biases()->get_grid_bias()},
        m_max_D_bias {params.m_max_D_bias} {

    // Setup grid bias (UGLY MESS)
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

    // Read in weights if specified
    if (params.m_biases_file != "") {
        read_weights(params.m_biases_file);
    }

    // Update starting configs if restarting
    if (m_params.m_restart_traj_file != "") {
        set_config_from_traj(m_params.m_restart_traj_file, params.m_restart_step);
    }
}

USGCMCSimulation::~USGCMCSimulation() {
}

void USGCMCSimulation::run() {
    int n {0};
    run_equilibration();
    while (n != m_max_num_iters) {
        run_iteration(n);
        n++;
        if (weights_converged()) {
            break;
        }
    }
    run_production(n);
}

void USGCMCSimulation::run_equilibration() {
    string postfix {"_iter-equil"};
    string output_filebase {m_params.m_output_filebase + postfix};
    m_output_files = setup_output_files(m_params, output_filebase, m_origami_system);
    m_logging_stream = new ofstream {output_filebase + ".out"};
    simulate(m_equil_steps);
    close_output_files();
    delete m_logging_stream;
}

void USGCMCSimulation::run_iteration(int n) {

    // Write each iteration's output to a seperate file
    string prefix {"_iter-" + std::to_string(n)};
    string output_filebase {m_params.m_output_filebase + prefix};

    bool step_is_equil {true};
    while (step_is_equil) {
        m_output_files = setup_output_files(m_params, output_filebase, m_origami_system);
        m_logging_stream = new ofstream {output_filebase + ".out"};

        clear_grids();
        simulate(m_steps);
        fill_grid_sets();
        step_is_equil = iteration_equilibrium_step();
        if (step_is_equil) {
            *m_us_stream << "Discarding iteration\n\n";
            m_s_i.clear();
            m_f_i.clear();
        }
        else {
            m_S_n.insert(m_s_i.begin(), m_s_i.end());
            estimate_current_weights();
            update_grids(n);
            update_bias(n);
            output_summary(n);
            close_output_files();
            delete m_logging_stream;
        }
    }
    output_weights();
}

void USGCMCSimulation::clear_grids() {
    m_points.clear();
    m_s_i.clear();
    m_f_i.clear();
    m_new_points.clear();
    m_old_points.clear();
    m_old_only_points.clear();
}

bool USGCMCSimulation::weights_converged() {

    // For now not using this
    return false;

    /*bool converged {true};
    if (m_old_lP_n.size() != m_S_n.size()) {
        converged = false;
        return converged;
    }

    // Check states that were sampled previously are sampled now
    for (auto point: m_S_n) {
        if (m_w_i[point] == 0) {
            converged = false;
            return converged;
        }
    }

    // Check that relative change in weight is less than max value
    // Maybe only check weights that aren't very small
    for (auto point: m_S_n) {
        double rel_dif {(m_lP_n[point] - m_old_lP_n[point]) / m_old_lP_n[point]};
        if (abs(rel_dif) > m_max_rel_P_diff) {
            converged = false;
                break;
        }
    }

    return converged;
    */
}

void USGCMCSimulation::run_production(int n) {
    // Hacky way to get relative weights right
    m_steps = m_prod_steps;

    string postfix {"_iter-prod"};
    string output_filebase {m_params.m_output_filebase + postfix};
    m_output_files = setup_output_files(m_params, output_filebase, m_origami_system);
    m_logging_stream = new ofstream {output_filebase + ".out"};

    simulate(m_prod_steps);
    fill_grid_sets();
    m_S_n.insert(m_s_i.begin(), m_s_i.end());
    estimate_current_weights();
    update_grids(n);
    prod_output_summary();
    close_output_files();
    delete m_logging_stream;
}

vector<GridPoint> USGCMCSimulation::get_points() {
    return m_points;
}

void USGCMCSimulation::read_weights(string filename) {
    ifstream jsonraw {filename, ifstream::binary};
    Json::Value jsonroot;
    jsonraw >> jsonroot;

    for (auto json_entry:jsonroot["biases"]) {
        GridPoint point {};
        for (auto comp: json_entry["point"]) {
            point.push_back(comp.asInt());
        }
        double point_bias {json_entry["bias"].asDouble()};
        m_E_w[point] = point_bias;
    }
}

void USGCMCSimulation::output_weights() {
    string filename {m_params.m_output_filebase + ".biases"};
    ofstream weights_file {filename};
    weights_file << "{\n";
    weights_file << "    \"biases\": [\n";
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
    weights_file << "    ]\n";
    weights_file << "}\n";
}

void USGCMCSimulation::set_config_from_traj(string filename, int step) {
    OrigamiTrajInputFile traj_inp {filename};
    Chains config {traj_inp.read_config(step)};
    m_origami_system.set_config(config);
}

void USGCMCSimulation::set_output_stream(ostream* out_stream) {
    m_us_stream = out_stream;
}

int USGCMCSimulation::get_grid_dim() {
    return m_grid_params.size();
}

void USGCMCSimulation::update_internal(long int step) {
    GridPoint point {};
    for (auto grid_param: m_grid_params) {
        point.push_back(grid_param->get_param());
    }
    m_s_i.insert(point);
    m_f_i[point] ++;

    int op_freq {m_params.m_order_params_output_freq};
    if (op_freq != 0 and step % op_freq == 0) {
        m_points.push_back(point);
    }
}

void USGCMCSimulation::estimate_current_weights() {
    // Average bias weights of iteration
    double ave_bias_weight {0};
    for (auto point: m_s_i) {
        double point_bias {m_grid_bias->calc_bias(point)};
        ave_bias_weight += m_f_i[point] * std::exp(point_bias);
    }

    // Calculate for all visited points
    for (auto point: m_s_i) {
        double point_bias {m_grid_bias->calc_bias(point)};
        double p_n_k {m_f_i[point] * std::exp(point_bias) /
            ave_bias_weight};
        m_p_i[point] = p_n_k;
    }

    // Update vector for unvisted points
    for (auto point: m_old_only_points) {
        m_p_i[point] = 0;
    }
}

bool SimpleUSGCMCSimulation::iteration_equilibrium_step() {
    return false;
}

void USGCMCSimulation::fill_grid_sets() {

    m_new_points.resize(m_s_i.size());
    m_old_points.resize(m_s_i.size());
    m_old_only_points.resize(m_S_n.size());

    // First find states that were visited before and now
    vector<GridPoint>::iterator it;
    it = std::set_intersection(m_s_i.begin(), m_s_i.end(),
            m_S_n.begin(), m_S_n.end(), m_old_points.begin());
    m_old_points.resize(it - m_old_points.begin());

    // Now find newly sampled points
    it = std::set_difference(m_s_i.begin(), m_s_i.end(), m_old_points.begin(),
            m_old_points.end(), m_new_points.begin());
    m_new_points.resize(it - m_new_points.begin());

    // Now find points sampled before not in current
    it = std::set_difference(m_S_n.begin(), m_S_n.end(), m_old_points.begin(),
            m_old_points.end(), m_old_only_points.begin());
    m_old_only_points.resize(it - m_old_only_points.begin());

    return;
}

void USGCMCSimulation::prod_output_summary() {
    *m_us_stream << "Production run\n";
    *m_us_stream << "Gridpoint w, p, P, E:\n";
    for (auto point: m_S_n) {
        for (auto coor: point) {
            *m_us_stream << coor << " ";
        }
        *m_us_stream << std::setprecision(3);
        *m_us_stream << ": " << std::setw(10) << m_w_i[point] << std::setw(10) <<
                m_p_i[point] << std::setw(10) << m_lP_n[point] <<
                std::setw(10) << m_E_w[point] << "\n";
    }
    *m_us_stream << "\n";
}

SimpleUSGCMCSimulation::SimpleUSGCMCSimulation(OrigamiSystem& origami,
        InputParameters& params):
        USGCMCSimulation(origami, params) {
}

void SimpleUSGCMCSimulation::update_bias(int) {
    m_old_lP_n = m_lP_n;
    for (auto point: m_S_n) {

        // Update big P
        m_lP_n[point] = m_p_i[point];

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
    m_grid_bias->replace_biases(m_E_w);
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

MWUSGCMCSimulation::MWUSGCMCSimulation(OrigamiSystem& origami,
        InputParameters& params) :
        GCMCSimulation {origami, params},
        m_params {params},
        m_max_num_iters {params.m_max_num_iters},
        m_system_order_params {dynamic_cast<OrigamiSystemWithBias*>(&origami)->
                get_system_order_params()},
        m_system_biases {dynamic_cast<OrigamiSystemWithBias*>(&origami)->
                get_system_biases()} {

    parse_windows_file(params.m_windows_file);

    // Prepare master node variables
    for (int i {0}; i != m_windows; i++) {
        m_points.push_back({});
        m_sims_converged.push_back(false);
        m_starting_files.push_back("");
        m_starting_steps.push_back(0);
        m_current_iters.push_back(0);

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
        string output_filebase {params.m_output_filebase + window_postfix};
        m_output_filebases.push_back(output_filebase);
    }

    // Setup square well potentials
    GridPoint window_min {m_window_mins[m_rank]};
    GridPoint window_max {m_window_maxs[m_rank]};
    OrderParam* n_staples {&m_system_order_params->get_num_staples()};
    OrderParam* n_domains {&m_system_order_params->get_num_bound_domains()};

    // Hard coded which ops are used
    int domain_min {window_min[0]};
    int domain_max {window_max[0]};
    m_system_biases->add_square_well_bias(n_domains, domain_min, domain_max,
                params.m_well_bias, params.m_outside_bias);
    int staple_min {window_min[1]};
    int staple_max {window_max[1]};
    m_system_biases->add_square_well_bias(n_staples, staple_min, staple_max,
                params.m_well_bias, params.m_outside_bias);

    // Update filebases and construct US sim object
    string window_postfix {m_window_postfixes[m_rank]};
    string output_filebase {m_output_filebases[m_rank]};
    params.m_output_filebase = output_filebase;
    if (params.m_biases_filebase != "") {
        params.m_biases_file = params.m_biases_filebase + window_postfix;
        params.m_biases_file += ".biases";
    }
    if (params.m_restart_traj_filebase != "") {
        params.m_restart_traj_filebase += window_postfix;
        params.m_restart_traj_file = params.m_restart_traj_filebase +
                m_params.m_restart_traj_postfix;
    }
    if (not params.m_restart_traj_files.empty()) {
        params.m_restart_traj_file = params.m_restart_traj_files[m_rank];
        params.m_restart_step = params.m_restart_steps[m_rank];
    }

    m_us_sim = new SimpleUSGCMCSimulation {origami, params};
    m_us_stream = new ofstream {output_filebase + ".out"};
    m_us_sim->set_output_stream(m_us_stream);
    m_grid_dim = m_us_sim->get_grid_dim();

}

MWUSGCMCSimulation::~MWUSGCMCSimulation() {
    delete m_us_sim;
    delete m_us_stream;
}

void MWUSGCMCSimulation::run() {
    m_us_sim->run_equilibration();
    bool sim_converged {false};
    int n {0};
    int n_sims {0};

    // This is really ugly
    while (n != m_max_num_iters) {
        if (m_rank != m_master_node and sim_converged) {
            break;
        }
        else if (m_rank == m_master_node and 
                std::all_of(m_sims_converged.begin(), m_sims_converged.end(),
                        [](bool i){return i;})) {
            break;
        }
        if (not sim_converged) {
            m_us_sim->run_iteration(n);
            sim_converged = m_us_sim->weights_converged();
        }
        update_master_order_params(n);
        update_master_converged_sims(sim_converged, n);
        update_starting_config(n);
        n++;
        if (not sim_converged) {
            n_sims++;
        }
        if (m_rank == m_master_node) {
            for (int i {0}; i != m_windows; i++) {
                if (not m_sims_converged[i]) {
                    m_current_iters[i]++;
                }
            }
            cout << "Iteration: " << n << "\n";
            cout << "Simulations converged: " <<  m_num_sims_converged << "/" <<
                m_windows << "\n";
        }
    }

    n_sims++;
    m_us_sim->run_production(n_sims);
}

void MWUSGCMCSimulation::update_master_order_params(int n) {
    if (m_rank != m_master_node) {

        // Must flatten to send
        vector<int> flattened_points {};
        vector<GridPoint> points {m_us_sim->get_points()};
        for (auto point: points) {
            for (auto comp: point) {
                flattened_points.push_back(comp);
            }
        }
        m_world.send(m_master_node, n, flattened_points);
    }
    else {
        m_points[m_master_node] = m_us_sim->get_points();
        for (int i {1}; i != m_windows; i++) {
            if (m_sims_converged[i]) {
                continue;
            }
            vector<int> f_points;
            m_world.recv(i, n, f_points);

            // Unflatten
            vector<GridPoint> points {};
            size_t j {0};
            while (j != f_points.size()) {
                points.push_back({});
                for (int k {0}; k != m_grid_dim; k++) {
                    points.back().push_back(f_points[j]);
                    j++;
                }
            }
            m_points[i] = points;
        }
    }
}

void MWUSGCMCSimulation::update_master_converged_sims(bool sim_converged, int n) {
    if (m_rank != m_master_node) {
        m_world.send(m_master_node, n, sim_converged);
    }
    else {
        if (not m_sims_converged[m_master_node]) {
            m_sims_converged[m_master_node] = sim_converged;
            m_num_sims_converged += sim_converged;
        }
        for (int i {1}; i != m_windows; i++) {
            if (m_sims_converged[i]) {
                continue;
            }
            bool slave_sim_converged;
            m_world.recv(i, n, slave_sim_converged);
            m_sims_converged[i] = slave_sim_converged;
            m_num_sims_converged += slave_sim_converged;
        }
    }
}

void MWUSGCMCSimulation::update_starting_config(int n) {
    if (m_rank == m_master_node) {
        select_starting_configs();
        m_starting_file = m_starting_files[m_master_node];
        m_starting_step = m_starting_steps[m_master_node];
        for (int i {1}; i != m_windows; i++) {
            m_world.send(i, n, m_starting_files[i]);
            m_world.send(i, n, m_starting_steps[i]);
        }
        m_us_sim->set_config_from_traj(m_starting_file, m_starting_step);
    }
    else {
        m_world.recv(m_master_node, n, m_starting_file);
        m_world.recv(m_master_node, n, m_starting_step);
        m_us_sim->set_config_from_traj(m_starting_file, m_starting_step);
    }
}

void MWUSGCMCSimulation::select_starting_configs() {
    m_order_param_to_configs.clear();
    for (int i {0}; i != m_windows; i++) {
        for (size_t step {0}; step != m_points[i].size(); step++) {
            GridPoint point {m_points[i][step]};

            // Make this more robust
            GridPoint win_point {point[0], point[1]};
            if (m_order_param_to_configs.find(win_point) != m_order_param_to_configs.end()) {
                m_order_param_to_configs[win_point].push_back({i, step});
            }
            else {
                m_order_param_to_configs[win_point] = {{i, step}};
            }
        }
    }

    for (int i {0}; i != m_windows; i++) {
        if (m_sims_converged[i]) {
            continue;
        }
        GridPoint point {};
        while (m_order_param_to_configs.find(point) == m_order_param_to_configs.end()) {
            point.clear();
            for (size_t j {0}; j != m_window_mins[0].size(); j++) {
                int min_comp {m_window_mins[i][j]};
                int max_comp {m_window_maxs[i][j]};
                point.push_back(m_random_gens.uniform_int(min_comp, max_comp));
            }
        }
        vector<pair<int, int>> possible_configs {m_order_param_to_configs[point]};
        int sel_i {m_random_gens.uniform_int(0, possible_configs.size() - 1)};
        pair<int, int> selected_config {possible_configs[sel_i]};
        int window {selected_config.first};
        int current_iter {m_current_iters[window]};
        string filename {m_output_filebases[selected_config.first] + "_iter-" +
                std::to_string(current_iter) + ".trj"};
        m_starting_files[i] = filename;
        m_starting_steps[i] = selected_config.second;
    }
}

void MWUSGCMCSimulation::parse_windows_file(string filename) {
    ifstream file {filename};
    string window_raw;
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
