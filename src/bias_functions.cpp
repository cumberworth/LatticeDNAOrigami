// bias_functions.cpp

#include "bias_functions.h"
#include "files.h"

namespace biasFunctions {

    using std::cout;
    using std::set;

    using files::OrigamiBiasFunctionsFile;
    using orderParams::DistOrderParam;

    double BiasFunction::get_bias() {
        return m_bias;
    }

    LinearFunction::LinearFunction(OrderParam& order_param, double slope):
            m_order_param {order_param},
            m_slope {slope} {

        m_slope = slope;
        update_bias();
    }

    double LinearFunction::update_bias() {
        double new_bias;
        int param {m_order_param.get_param()};
        if (m_order_param.defined()) {
            new_bias = calc_bias(param);
        }
        else {
            new_bias = 0;
        }
        m_bias = new_bias;
        
        return new_bias;
    }

    double LinearFunction::calc_bias(int param) {
        double new_bias {0};
        new_bias = m_slope*param;
        return new_bias;
    }

    double LinearFunction::check_bias() {
        double checked_bias;
        int param {m_order_param.get_checked_param()};
        if (m_order_param.defined()) {
            checked_bias = calc_bias(param);
        }
        else {
            checked_bias = 0;
        }
        return checked_bias;
    }

    LinearStepBiasFunction::LinearStepBiasFunction(OrderParam& order_param,
            int min_param, int max_param, double max_bias) :
            m_order_param {order_param},
            m_min_param {min_param},
            m_max_param {max_param},
            m_max_bias {max_bias} {

        m_slope = m_max_bias / (max_param - min_param);
        update_bias();
    }

    double LinearStepBiasFunction::update_bias() {
        double new_bias;
        int param {m_order_param.get_param()};
        if (m_order_param.defined()) {
            new_bias = calc_bias(param);
        }
        else {
            new_bias = 0;
        }
        m_bias = new_bias;
        
        return new_bias;
    }

    double LinearStepBiasFunction::calc_bias(int param) {
        double new_bias {0};
        if (param <= m_min_param) {
            new_bias = 0;
        }
        else if (param > m_min_param and param <= m_max_param) {
            new_bias = m_slope * (param - m_min_param);
        }
        else {
            new_bias = m_max_bias;
        }
        return new_bias;
    }

    double LinearStepBiasFunction::check_bias() {
        double checked_bias;
        int param {m_order_param.get_checked_param()};
        if (m_order_param.defined()) {
            checked_bias = calc_bias(param);
        }
        else {
            checked_bias = 0;
        }
        return checked_bias;
    }

    LinearStepWellBiasFunction::LinearStepWellBiasFunction(
            OrderParam& order_param, int min_param, int max_param,
            double well_bias, double min_bias, double slope) :
            m_order_param {order_param},
            m_min_param {min_param},
            m_max_param {max_param},
            m_well_bias {well_bias},
            m_min_bias {min_bias},
            m_slope {slope} {

        update_bias();
    }

    double LinearStepWellBiasFunction::update_bias() {
        double new_bias;
        int param {m_order_param.get_param()};
        if (m_order_param.defined()) {
            new_bias = calc_bias(param);
        }
        else {
            new_bias = 0;
        }
        m_bias = new_bias;
        
        return new_bias;
    }

    double LinearStepWellBiasFunction::calc_bias(int param) {
        double new_bias {0};
        if (param < m_min_param) {
            new_bias = m_slope * (m_min_param - param - 1) + m_min_bias;
        }
        else if (param > m_max_param) {
            new_bias = m_slope * (param - m_max_param - 1) + m_min_bias;
        }
        else {
            new_bias = m_well_bias;
        }
        return new_bias;
    }

    double LinearStepWellBiasFunction::check_bias() {
        double checked_bias;
        int param {m_order_param.get_checked_param()};
        if (m_order_param.defined()) {
            checked_bias = calc_bias(param);
        }
        else {
            checked_bias = 0;
        }
        return checked_bias;
    }

    SquareWellBiasFunction::SquareWellBiasFunction(OrderParam& order_param,
            int min_param, int max_param, double well_bias, double outside_bias) :
            m_order_param {order_param},
            m_min_param {min_param},
            m_max_param {max_param},
            m_well_bias {well_bias},
            m_outside_bias {outside_bias} {

        update_bias();
    }

    double SquareWellBiasFunction::update_bias() {
        double new_bias;
        int param {m_order_param.get_param()};
        if (m_order_param.defined()) {
            new_bias = calc_bias(param);
        }
        else {
            new_bias = 0;
        }
        m_bias = new_bias;
        
        return new_bias;
    }

    double SquareWellBiasFunction::calc_bias(int param) {
        double new_bias {0};
        if (param < m_min_param or param > m_max_param) {
            new_bias = m_outside_bias;
        }
        else {
            new_bias = m_well_bias;
        }

        return new_bias;
    }

    double SquareWellBiasFunction::check_bias() {
        double checked_bias;
        int param {m_order_param.get_checked_param()};
        if (m_order_param.defined()) {
            checked_bias = calc_bias(param);
        }
        else {
            checked_bias = 0;
        }
        return checked_bias;
    }

    void SquareWellBiasFunction::set_min_op(int min_op) {
        m_min_param = min_op;
    }

    void SquareWellBiasFunction::set_max_op(int max_op) {
        m_max_param = max_op;
    }

    GridBiasFunction::GridBiasFunction(
            vector<reference_wrapper<OrderParam>> ops):

            m_ops {ops} {
    }

    int GridBiasFunction::get_dim() {
        return m_ops.size();
    }

    vector<int> GridBiasFunction::get_point() {
        vector<int> point {};
        for (auto op: m_ops) {
            point.push_back(op.get().get_param());
        }

        return point;
    }

    void GridBiasFunction::replace_biases(unordered_map<vector<int>, double> bias_grid) {
        m_bias_grid = bias_grid;
    }

    double GridBiasFunction::calc_bias(vector<int> params) {
        auto key_value = m_bias_grid.find(params);
        double g_bias;
        if (key_value == m_bias_grid.end()) {
            g_bias = m_off_grid_bias;
        }
        else {
            g_bias = key_value->second;
        }
        return g_bias;
    }

    double GridBiasFunction::update_bias() {
        vector<int> key {};
        for (auto op: m_ops) {

            // No bias if order parameters not all defined
            if (not op.get().defined()) {
                m_bias = 0;
                return 0;
            }
            key.push_back(op.get().get_param());
        }
        double bias {calc_bias(key)};
        m_bias = bias;

        return bias;
    }

    double GridBiasFunction::check_bias() {
        vector<int> key {};
        for (auto op: m_ops) {

            // No bias if order parameters not all defined
            if (not op.get().defined()) {
                return 0;
            }
            key.push_back(op.get().get_checked_param());
        }
        double checked_bias {calc_bias(key)};

        return checked_bias;
    }

    SystemBiases::SystemBiases(
            OrigamiSystem& origami,
            SystemOrderParams& system_order_params,
            InputParameters& params) :
            m_origami {origami},
            m_system_order_params {system_order_params},
            m_bias_mult {params.m_bias_funcs_mult} {

        // Setup dependency table
        vector<pair<int, int>> keys {};
        for (auto chain: origami.get_chains()) {
            for (auto domain: chain) {
                pair<int, int> key {domain->m_c, domain->m_d};
                m_domain_update_biases[key] = {};
                keys.push_back(key);
            }
        }

        if (params.m_bias_funcs_filename != "") {
            setup_biases(params.m_bias_funcs_filename, keys);
        }
    }

    void SystemBiases::setup_biases(
            string biases_filename,
            vector<pair<int, int>> keys) {

        OrigamiBiasFunctionsFile biases_file {biases_filename};
        vector<vector<string>> level_to_types {biases_file.get_types_by_level()};
        vector<vector<string>> level_to_tags {biases_file.get_tags_by_level()};
        vector<vector<string>> level_to_labels {biases_file.get_labels_by_level()};
        vector<vector<vector<string>>> level_to_ops {biases_file.get_ops_by_level()};
        vector<vector<vector<string>>> level_to_d_biases {biases_file.get_d_biases_by_level()};
        unordered_map<string, vector<pair<int, int>>> tag_to_domains {};
        for (size_t i {0}; i != level_to_types.size(); i++) {
            m_level_to_biases.push_back({});
            m_move_update_biases.push_back({});
            m_levels++;
            for (auto key: keys) {
                m_domain_update_biases[key].push_back({});
            }
            for (size_t j {0}; j!= level_to_types[i].size(); j++) {
                BiasFunction* bias_f;
                string type {level_to_types[i][j]};
                string tag {level_to_tags[i][j]};
                vector<string> op_tags {level_to_ops[i][j]};
                vector<string> d_bias_tags {level_to_d_biases[i][j]};
                set<pair<int, int>> d_domains {};
                bool update_per_domain {true};
                for (auto op_tag: op_tags) {
                    vector<pair<int, int>> domains {m_system_order_params.
                        get_dependent_domains(op_tag)};
                    if (domains.size() == 0) {
                        update_per_domain = false;
                        break;
                    }
                    else {
                        for (auto d: domains) {
                            d_domains.insert(d);
                        }
                    }
                }
                for (auto d_bias_tag: d_bias_tags) {
                    vector<pair<int, int>> domains {tag_to_domains[d_bias_tag]};
                    if (domains.size() == 0 or not update_per_domain) {
                        update_per_domain = false;
                        break;
                    }
                    for (auto d: domains) {
                        d_domains.insert(d);
                    }
                }
                tag_to_domains[tag] = {};
                if (update_per_domain) {
                    for (auto domain: d_domains) {
                        tag_to_domains[tag].push_back(domain);
                    }
                }
                if (type == "LinearStep") {
                    OrderParam& op {m_system_order_params.get_order_param(op_tags[0])};
                    int min_op {biases_file.get_int_option(i, j, "min_op")};
                    int max_op {biases_file.get_int_option(i, j, "max_op")};
                    double max_bias {biases_file.get_double_option(i, j, "max_bias")};
                    bias_f = new LinearStepBiasFunction {op, min_op,
                        max_op, max_bias};
                }
                if (type == "LinearStepWell") {
                    OrderParam& op {m_system_order_params.get_order_param(op_tags[0])};
                    int min_op {biases_file.get_int_option(i, j, "min_op")};
                    int max_op {biases_file.get_int_option(i, j, "max_op")};
                    double well_bias {biases_file.get_double_option(i, j, "well_bias")};
                    double min_bias {biases_file.get_double_option(i, j, "min_bias")};
                    double slope {biases_file.get_double_option(i, j, "slope")};
                    bias_f = new LinearStepWellBiasFunction {op, min_op,
                        max_op, well_bias, min_bias, slope};
                }
                else if (type == "SquareWell") {
                    OrderParam& op {m_system_order_params.get_order_param(op_tags[0])};
                    int min_op {biases_file.get_int_option(i, j, "min_op")};
                    int max_op {biases_file.get_int_option(i, j, "max_op")};
                    double well_bias {biases_file.get_double_option(i, j, "well_bias")};
                    double outside_bias {biases_file.get_double_option(i, j, "outside_bias")};
                    bias_f = new SquareWellBiasFunction {op,
                            min_op, max_op, well_bias, outside_bias};
                }
                else if (type == "Grid") {
                    vector<reference_wrapper<OrderParam>> ops {};
                    for (auto op_tag: op_tags) {
                        reference_wrapper<OrderParam> op {m_system_order_params.
                                get_order_param(op_tag)};
                        ops.push_back(op);
                    }
                    bias_f = new GridBiasFunction {ops};                        
                }
                else {
                    cout << "No such bias function type";
                    throw utility::SimulationMisuse {};
                }

                m_level_to_biases[i].emplace_back(bias_f);
                BiasFunction& bias_f_r {*m_level_to_biases[i].back()};
                m_tag_to_biases.emplace(tag, bias_f_r);
                if (update_per_domain) {
                    for (auto domain: d_domains) {
                        m_domain_update_biases.at(domain)[i].emplace_back(bias_f_r);
                    }
                    m_domain_update_bias += bias_f_r.get_bias();
                }
                else {
                    m_move_update_biases[i].emplace_back(bias_f_r);
                    m_move_update_bias += bias_f_r.get_bias();
                }
            }
        }
    }

    //void SystemBiases::add_chain(vector<Domain*> chain) {
        // For each domain
        // If biases that include all chains of type x
            // Create new order parameters and bias functions
            // Add dependencies
    //}

    //void SystemBiases::remove_chain(vector<Domain*> chain) {
        // For each domain
        // Find order params and bias fs dependent on domain
        // Remove dependecies
        // Iterate through order params and bias fs
        //     delete those with 0 dependencies
    //}

    double SystemBiases::calc_total() {
        double total_bias {0};
        for (auto& level: m_level_to_biases) {
            for (auto &bias_f: level) {
                double bias {bias_f->update_bias()};
                total_bias += bias;
            }
        }

        return total_bias * m_bias_mult;
    }

    double SystemBiases::get_total_bias() {
        return (m_domain_update_bias + m_move_update_bias) * m_bias_mult;
    }

    double SystemBiases::get_domain_update_bias() {
        return m_domain_update_bias * m_bias_mult;
    }

    double SystemBiases::get_move_update_bias() {
        return m_move_update_bias * m_bias_mult;
    }

    void SystemBiases::update_bias_mult(double bias_mult) {
        m_bias_mult = bias_mult;
    }

    double SystemBiases::calc_move() {
        double bias_diff {0};
        for (auto level: m_move_update_biases) {
            for (auto bias_f: level) {
                double prev_bias {bias_f.get().get_bias()};
                double new_bias {bias_f.get().update_bias()};
                bias_diff += new_bias - prev_bias;
            }
        }
        m_move_update_bias += bias_diff;

        return bias_diff * m_bias_mult;
    }

    double SystemBiases::calc_one_domain(Domain& domain) {
        double bias_diff {0};

        // Get dependent bias functions and update
        pair<int, int> key {domain.m_c, domain.m_d};
        for (auto level: m_domain_update_biases[key]) {
            for (auto bias_f: level) {
                double prev_bias {bias_f.get().get_bias()};
                double new_bias {bias_f.get().update_bias()};
                bias_diff += new_bias - prev_bias;
            }
        }
        m_domain_update_bias += bias_diff;

        return bias_diff * m_bias_mult;
    }

    double SystemBiases::check_one_domain(Domain& domain) {
        double bias_diff {0};

        // Get dependent bias functions
        pair<int, int> key {domain.m_c, domain.m_d};
        for (auto level: m_domain_update_biases[key]) {
            for (auto bias_f: level) {

                // Check bias bias values without internal update
                double prev_bias {bias_f.get().get_bias()};
                double new_bias {bias_f.get().check_bias()};
                bias_diff += new_bias - prev_bias;
            }
        }

        return bias_diff * m_bias_mult;
    }

    GridBiasFunction& SystemBiases::get_grid_bias(string tag) {
        return static_cast<GridBiasFunction&>(m_tag_to_biases.at(tag).get());
    }

    SquareWellBiasFunction& SystemBiases::get_square_well_bias(string tag) {
        return static_cast<SquareWellBiasFunction&>(m_tag_to_biases.at(tag).get());
    }
}
