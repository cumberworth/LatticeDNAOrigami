// us_simulation.cpp

#include <cmath>
#include <cstdlib>

#include "json/json.h"

#include "files.h"
#include "us_simulation.h"
#include "utility.h"

namespace us {

    using std::ifstream;

    using biasFunctions::SquareWellBiasFunction;
    using files::OrigamiTrajInputFile;
    using origami::Chains;
    using utility::SimulationMisuse;

    USGCMCSimulation::USGCMCSimulation(
            OrigamiSystem& origami,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params) :
            GCMCSimulation(origami, ops, biases, params),
            m_params {params},
            m_max_num_iters {params.m_max_num_iters},
            m_equil_steps {params.m_equil_steps},
            m_max_equil_dur {params.m_max_equil_dur},
            m_iter_steps {params.m_iter_steps},
            m_max_iter_dur {params.m_max_iter_dur},
            m_prod_steps {params.m_prod_steps},
            m_max_prod_dur {params.m_max_prod_dur},
            m_max_rel_P_diff {params.m_max_rel_P_diff},
            m_grid_bias {biases.get_grid_bias(params.m_us_grid_bias_tag)},
            m_max_D_bias {params.m_max_D_bias} {

        // Read in weights if specified
        if (params.m_biases_file != "") {
            ifstream bias_file {params.m_biases_file};
            if (bias_file.good()) {
                read_weights(params.m_biases_file);
            }
        }

        // Update starting configs if restarting
        if (m_params.m_restart_traj_file != "") {
            set_config_from_traj(m_params.m_restart_traj_file, params.m_restart_step);
        }
    }

    void USGCMCSimulation::run() {
        int n;
        run_equilibration();
        for (n = 0; n != m_max_num_iters; n++) {
            run_iteration(n);
        }
        run_production(n);
    }

    void USGCMCSimulation::run_equilibration() {
        set_max_dur(m_max_equil_dur);

        // Setup output files
        string postfix {"_iter-equil"};
        string output_filebase {m_params.m_output_filebase + postfix};
        m_output_files = simulation::setup_output_files(m_params,
                output_filebase, m_origami_system, m_ops, m_biases);
        m_logging_stream = new ofstream {output_filebase + ".out"};

        m_steps = m_equil_steps;
        m_steps = simulate(m_steps);

        // Cleanup
        close_output_files();
        delete m_logging_stream;
    }

    void USGCMCSimulation::run_iteration(int n) {
        set_max_dur(m_max_iter_dur);

        // Write each iteration's output to a seperate file
        string prefix {"_iter-" + std::to_string(n)};
        string output_filebase {m_params.m_output_filebase + prefix};
        m_output_files = simulation::setup_output_files(m_params,
                output_filebase, m_origami_system, m_ops, m_biases);
        m_logging_stream = new ofstream {output_filebase + ".out"};

        m_steps = m_iter_steps;
        clear_grids();
        m_steps = simulate(m_steps);
        fill_grid_sets();
        m_S_n.insert(m_s_i.begin(), m_s_i.end());
        estimate_current_weights();
        update_grids(n);
        update_bias(n);
        output_summary(n);

        // Cleanup
        close_output_files();
        delete m_logging_stream;
        output_weights();
    }

    void USGCMCSimulation::clear_grids() {
        m_points.clear();
        m_s_i.clear();
        m_f_i.clear();
        m_w_i.clear();
        m_new_points.clear();
        m_old_points.clear();
        m_old_only_points.clear();
    }

    void USGCMCSimulation::run_production(int n) {
        set_max_dur(m_max_prod_dur);

        // Setup output files
        string postfix {"_iter-prod"};
        string output_filebase {m_params.m_output_filebase + postfix};
        m_output_files = simulation::setup_output_files(m_params,
                output_filebase, m_origami_system, m_ops, m_biases);
        m_logging_stream = new ofstream {output_filebase + ".out"};

        m_steps = m_prod_steps;
        clear_grids();
        m_steps = simulate(m_steps);
        fill_grid_sets();
        m_S_n.insert(m_s_i.begin(), m_s_i.end());
        estimate_current_weights();
        update_grids(n);
        output_summary(n);

        // Cleanup
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

    void USGCMCSimulation::set_output_stream(ostream* out_stream) {
        m_us_stream = out_stream;
    }

    int USGCMCSimulation::get_grid_dim() {
        return m_grid_bias.get_dim();
    }

    void USGCMCSimulation::update_internal(long long int step) {
        GridPoint point {m_grid_bias.get_point()};
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
            double point_bias {m_grid_bias.calc_bias(point)};
            ave_bias_weight += m_f_i[point] * std::exp(point_bias);
        }

        // Calculate for all visited points
        for (auto point: m_s_i) {
            double point_bias {m_grid_bias.calc_bias(point)};
            double p_n_k {m_f_i[point] * std::exp(point_bias) /
                ave_bias_weight};
            m_p_i[point] = p_n_k;
        }

        // Update vector for unvisted points
        for (auto point: m_old_only_points) {
            m_p_i[point] = 0;
        }
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

    SimpleUSGCMCSimulation::SimpleUSGCMCSimulation(
            OrigamiSystem& origami,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params):
            USGCMCSimulation(origami, ops, biases, params) {
    }

    void SimpleUSGCMCSimulation::update_bias(int) {
        m_old_p_i = m_p_i;
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
            SystemBiases & biases,
            InputParameters& params) :
            GCMCSimulation {origami, ops, biases, params},
            m_params {params},
            m_max_num_iters {params.m_max_num_iters},
            m_local_dir {params.m_local_dir},
            m_central_dir {params.m_central_dir} {

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
            m_us_sim->run_iteration(n);
            copy_files_to_central_dir(n);
            update_master_order_params(n);
            update_starting_config(n);
            if (m_rank == m_master_node) {
                output_iter_summary(n);
            }
        }

        m_us_sim->run_production(n);
    }

    void MWUSGCMCSimulation::setup_window_variables() {
        for (int i {0}; i != m_windows; i++) {
            m_points.push_back({});
            m_starting_files.push_back("");
            m_starting_steps.push_back(0);

            // Filename bases
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
            string output_filebase {m_params.m_output_filebase +
                    window_postfix};
            m_output_filebases.push_back(output_filebase);

            // Calculate number of grid points in the window
            unsigned int num_points {1};
            for (size_t j {0}; j != m_window_mins[i].size(); j++) {
                num_points *= (m_window_maxs[i][j] - m_window_mins[i][j] + 1);
            }
            m_num_points.push_back(num_points);
        }
    }

    void MWUSGCMCSimulation::setup_window_restraints() {
        GridPoint window_min {m_window_mins[m_rank]};
        GridPoint window_max {m_window_maxs[m_rank]};

        for (size_t i {0}; i != m_window_bias_tags.size(); i++) {
            string bias_tag {m_window_bias_tags[i]};
            SquareWellBiasFunction& bias_f {m_biases.get_square_well_bias(
                    bias_tag)};
            int op_min {window_min[i]};
            int op_max {window_max[i]};
            bias_f.set_min_op(op_min);
            bias_f.set_max_op(op_max);
        }
    }

    void MWUSGCMCSimulation::setup_window_sims(OrigamiSystem& origami) {

        // US sims reads filebase names from params object, so update
        string window_postfix {m_window_postfixes[m_rank]};
        string output_filebase {m_local_dir + "/" + m_output_filebases[m_rank]};
        m_params.m_output_filebase = output_filebase;

        // If available read modify filename for input biases
        if (m_params.m_biases_filebase != "") {
            m_params.m_biases_file = m_params.m_biases_filebase +
                    window_postfix;
            m_params.m_biases_file += ".biases";
        }

        // If available modify filename for seperate starting config for each window
        // Standard names
        if (m_params.m_restart_traj_filebase != "") {
            m_params.m_restart_traj_filebase += window_postfix;
            m_params.m_restart_traj_file = m_params.m_restart_traj_filebase +
                    m_params.m_restart_traj_postfix;
        }

        // Names invididually specified in param file
        if (not m_params.m_restart_traj_files.empty()) {
            m_params.m_restart_traj_file =
                    m_params.m_restart_traj_files[m_rank];
            m_params.m_restart_step = m_params.m_restart_steps[m_rank];
        }

        // Create simulation objects
        m_us_sim = new SimpleUSGCMCSimulation {origami, m_ops, m_biases,
                m_params};
        m_us_stream = new ofstream {output_filebase + ".out"};
        m_us_sim->set_output_stream(m_us_stream);
        m_grid_dim = m_us_sim->get_grid_dim();
    }

    void MWUSGCMCSimulation::output_iter_summary(int n) {
        cout << "Iteration: " << n << "\n";
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

    void MWUSGCMCSimulation::copy_files_to_central_dir(int n) {

        // Move most recent iteration's files to central dir
        string filebase_cur {m_output_filebases[m_rank] +
                "_iter-" + std::to_string(n)};
        string move_command {"mv " + m_local_dir + "/" + filebase_cur +
            "* " + m_central_dir + "/"};
        std::system(move_command.c_str());

        // Remove previous iteration files from central dir
        string filebase_prev {m_output_filebases[m_rank] +
                "_iter-" + std::to_string(n - 1)};
        string del_command {"rm -f " + m_central_dir + "/" + filebase_prev +
            "*"};
        std::system(del_command.c_str());
    }

    void MWUSGCMCSimulation::update_starting_config(int n) {
        if (m_rank == m_master_node) {
            select_starting_configs(n);
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

    void MWUSGCMCSimulation::select_starting_configs(int n) {
        sort_configs_by_ops();

        // For each window select an order parameter and then associated config
        for (int i {0}; i != m_windows; i++) {

            // Select a point in the window that has been sampled
            GridPoint point {};
            set<GridPoint> tried_points {};
            bool point_sampled {false};
            bool all_points_tried {false};
            while (not point_sampled) {
                point.clear();

                // Uniformly draw without replacement point in window's domain
                for (size_t j {0}; j != m_window_mins[0].size(); j++) {
                    int min_comp {m_window_mins[i][j]};
                    int max_comp {m_window_maxs[i][j]};
                    point.push_back(m_random_gens.uniform_int(min_comp,
                            max_comp));
                }
                point_sampled = (m_order_param_to_configs.find(point) !=
                        m_order_param_to_configs.end());
                if (point_sampled) {
                    break;
                }
                else {
                    tried_points.insert(point);
                    all_points_tried = (tried_points.size() == m_num_points[i]);
                    if (all_points_tried) {
                        cout << "Window " << i << " outside of domain\n";
                        throw SimulationMisuse {};
                    }
                }
            }

            // Uniformly select one of the configs with the selected op set
            vector<pair<int, int>> possible_configs {
                    m_order_param_to_configs[point]};
            int sel_i {m_random_gens.uniform_int(0, possible_configs.size() -
                    1)};
            pair<int, int> selected_config {possible_configs[sel_i]};
            string filename {m_central_dir + "/" +
                    m_output_filebases[selected_config.first] +
                    "_iter-" + std::to_string(n) + ".trj"};
            m_starting_files[i] = filename;
            m_starting_steps[i] = selected_config.second;
        }
    }

    void MWUSGCMCSimulation::sort_configs_by_ops() {

        // For each window order available configs by their order parameters
        m_order_param_to_configs.clear();
        for (int i {0}; i != m_windows; i++) {
            for (size_t step {0}; step != m_points[i].size(); step++) {

                // Would need to modify to allow the use of a subset of
                // the ops collected during the simulation
                GridPoint point {m_points[i][step]};
                bool point_sampled {m_order_param_to_configs.find(point) !=
                        m_order_param_to_configs.end()};
                if (point_sampled) {
                    m_order_param_to_configs[point].push_back({i, step});
                }
                else {
                    m_order_param_to_configs[point] = {{i, step}};
                }
            }
        }
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
}
