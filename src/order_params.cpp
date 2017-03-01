// order_params.cpp

#include "utility.h"
#include "order_params.h"

using namespace Utility;
using namespace OrderParams;

using std::abs;

DistOrderParam::DistOrderParam(Domain& domain_1, Domain& domain_2) :
        m_domain_1 {domain_1},
        m_domain_2 {domain_2} {
    calc_param();
}

int DistOrderParam::calc_param() {
    int dist {-1};
    if (m_domain_1.m_state != Occupancy::unassigned and
            m_domain_2.m_state != Occupancy::unbound) {

        m_defined = true;
        VectorThree diff_vec {m_domain_2.m_pos - m_domain_1.m_pos};
        dist = diff_vec.abssum();
    }
    else {
        m_defined = false;
    }
    m_param = dist;

    return dist;
}

int DistOrderParam::check_param(Domain& domain, VectorThree new_pos, VectorThree,
        Occupancy) {
    Domain* unmodded_domain {&m_domain_1};
    if (domain.m_d == m_domain_1.m_d) {
        unmodded_domain = &m_domain_2;
    }
    int dist {-1};
    if (unmodded_domain->m_state != Occupancy::unassigned) {
        dist = (new_pos - unmodded_domain->m_pos).abssum();
        // Consider using a checked defined variable instead
        m_defined = true;
    }
    else {
        m_defined = false;
    }

    return dist;
}

int DistOrderParam::get_param() {
    return m_param;
}

bool DistOrderParam::dependent_on(Domain& domain) {
    bool dependent {false};
    if (&domain == &m_domain_1 or &domain == &m_domain_2) {
        dependent = true;
    }

    return dependent;
}

vector<Domain*> DistOrderParam::get_depending_domains() {
    return {&m_domain_1, &m_domain_2};
}

bool DistOrderParam::defined() {
    return m_defined;
}

DistSumOrderParam::DistSumOrderParam(vector<OrderParam*> dist_params) :
        m_dist_params {dist_params} {
    calc_param();
}

int DistSumOrderParam::calc_param() {
    int dist_sum {0};
    for (auto param: m_dist_params) {
        dist_sum += param->get_param();
    }
    m_param = dist_sum;

    return m_param;
}

int DistSumOrderParam::get_param() {
    return m_param;
}

bool DistSumOrderParam::dependent_on(Domain& domain) {
    bool dependent {false};
    for (auto dist_param: m_dist_params) {
        dependent = dist_param->dependent_on(domain);
        if (dependent == true) {
            break;
        }
    }
    
    return dependent;
}

vector<Domain*> DistSumOrderParam::get_depending_domains() {
//    vector<Domain&> domains {};
//    for (auto dist_param: m_dist_params) {
//        dependent = dist_param.dependent_on(domain);
//        if (dependent == true) {
//            break;
//        }
    return {};
}

bool DistSumOrderParam::defined() {
    return true;
}

NumStaplesOrderParam::NumStaplesOrderParam(OrigamiSystem& origami) :
        m_origami {origami} {
}

int NumStaplesOrderParam::calc_param() {
    return m_origami.num_staples();
}

int NumStaplesOrderParam::get_param() {
    return m_origami.num_staples();
}

bool NumStaplesOrderParam::dependent_on(Domain& domain) {
    bool dependent;
    if (domain.m_c > 0) {
        dependent = true;
    }
    else {
        dependent = false;
    }

    return dependent;
}

vector<Domain*> NumStaplesOrderParam::get_depending_domains() {
    return {};
}

bool NumStaplesOrderParam::defined() {
    return true;
}

NumBoundDomainPairsOrderParam::NumBoundDomainPairsOrderParam(
        OrigamiSystem& origami) :
        m_origami {origami} {
}

int NumBoundDomainPairsOrderParam::calc_param() {
    return m_origami.num_staples();
}

int NumBoundDomainPairsOrderParam::get_param() {
    return m_origami.num_bound_domain_pairs();
}

bool NumBoundDomainPairsOrderParam::dependent_on(Domain& domain) {
    bool dependent;
    if (domain.m_c == 0) {
        dependent = true;
    }
    else {
        dependent = false;
    }

    return dependent;
}

vector<Domain*> NumBoundDomainPairsOrderParam::get_depending_domains() {
    return {};
}

bool NumBoundDomainPairsOrderParam::defined() {
    return true;
}

SystemOrderParams::SystemOrderParams(InputParameters& params,
        OrigamiSystem& origami) :
        m_origami {origami} {

    // Setup dependency table
    for (auto chain: origami.get_chains()) {
        for (auto domain: chain) {
            pair<int, int> key {domain->m_c, domain->m_d};
            m_domain_to_order_params[key] = {};
        }
    }

    // Setup each type of order parameter and bias function
    if (params.m_distance_pairs.size() != 0) {
        setup_distance_param(params);
    }
    if (params.m_distance_sum) {
        DistSumOrderParam* dist_sum;
        dist_sum = new DistSumOrderParam {get_distance_params()};
        m_order_params.push_back(dist_sum);
    }
}

vector<OrderParam*> SystemOrderParams::get_distance_params() {
    // For now only distance order params, but will need to change eventually
    return m_dist_order_params;
}

void SystemOrderParams::setup_distance_param(InputParameters& params) {

    // Extract domain ids and store domain pairs
    vector<pair<Domain&, Domain&>> domain_pairs {};
    size_t num_restrained_domains {params.m_distance_pairs.size()};
    for (size_t pair_i {0}; pair_i != num_restrained_domains; pair_i += 2) {
        int domain_1_i {params.m_distance_pairs[pair_i]};
        Domain& domain_1 {*m_origami.get_domain(0, domain_1_i)}; 
        int domain_2_i {params.m_distance_pairs[pair_i + 1]};
        Domain& domain_2 {*m_origami.get_domain(0, domain_2_i)};
        domain_pairs.push_back({domain_1, domain_2});
    }

    for (size_t pair_i {0}; pair_i != domain_pairs.size(); pair_i++) {
        Domain& domain_1 {domain_pairs[pair_i].first};
        Domain& domain_2 {domain_pairs[pair_i].second};
        OrderParam* order_param;
        order_param = new DistOrderParam {domain_1, domain_2};
        m_simple_order_params.push_back(order_param);
        add_param_dependency(domain_1, order_param);
        add_param_dependency(domain_2, order_param);
    }
}

void SystemOrderParams::add_param_dependency(Domain& domain, OrderParam* order_param) {
    int c_i {domain.m_c};
    int d_i {domain.m_d};
    pair<int, int> key {c_i, d_i};
    m_domain_to_order_params[key].push_back(order_param);
}

vector<OrderParam*> SystemOrderParams::get_dependent_order_params(Domain& domain) {
    pair<int, int> key {domain.m_c, domain.m_d};
    return m_domain_to_order_params[key];
}

void SystemOrderParams::update_one_domain(Domain& domain) {
    vector<OrderParam*> order_params {get_dependent_order_params(domain)};
    for (auto order_param: order_params) {
        order_param->calc_param();
    }
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

LinearStepBiasFunction::~LinearStepBiasFunction() {
    delete &m_order_param;
}

double LinearStepBiasFunction::update_bias() {
    double bias;
    int param {m_order_param.get_param()};
    if (m_order_param.defined()) {
        bias = calc_bias(param);
    }
    else {
        bias = 0;
    }
    m_bias = bias;
    
    return bias;
}

double LinearStepBiasFunction::calc_bias(int param) {
    double bias {0};
    if (param <= m_min_param) {
        bias = 0;
    }
    else if (param > m_min_param and param <= m_max_param) {
        bias = m_slope * (param - m_min_param);
    }
    else {
        bias = m_max_bias;
    }
    return bias;
}

double LinearStepBiasFunction::check_bias(Domain& domain, VectorThree new_pos,
        VectorThree new_ore, Occupancy state) {
    double bias;
    int param {m_order_param.check_param(domain, new_pos, new_ore, state)};
    if (m_order_param.defined()) {
        bias = calc_bias(param);
    }
    else {
        bias = 0;
    }
    return bias;
}

bool LinearStepBiasFunction::dependent_on(OrderParam& order_param) {
    bool dependent {false};
    if (&order_param == &m_order_param) {
        dependent = true;
    }

    return dependent;
}

double LinearStepBiasFunction::get_bias() {
    return m_bias;
}

BinBiasFunction::BinBiasFunction(OrderParam& order_param,
        unordered_map<int, double> bins_to_biases) :
        m_order_param {order_param},
        m_bins_to_biases {bins_to_biases} {
    update_bias();
}

double BinBiasFunction::calc_bias(int param) {
    return m_bins_to_biases[param];
}

double BinBiasFunction::update_bias() {
    int param {m_order_param.get_param()};
    double bias {calc_bias(param)};
    m_bias = bias;

    return bias;
}

bool BinBiasFunction::dependent_on(OrderParam& order_param) {
    bool dependent;
    if (&order_param == &m_order_param) {
        dependent = true;
    }
    else {
        dependent = false;
    }

    return dependent;
}

double BinBiasFunction::get_bias() {
    return m_bias;
}

SystemBiases::SystemBiases(OrigamiSystem& origami,
            SystemOrderParams& system_order_params, InputParameters& params) :
        m_origami {origami},
        m_system_order_params {system_order_params},
        m_bias_mult {params.m_bias_mult} {

    // Setup dependency table
    for (auto chain: origami.get_chains()) {
        for (auto domain: chain) {
            pair<int, int> key {domain->m_c, domain->m_d};
            m_domain_to_bias_fs[key] = {};
        }
    }

    // Setup each type of order parameter and bias function
    if (params.m_distance_bias) {
        setup_distance_bias(params);
    }
}

SystemBiases::~SystemBiases() {
    for (auto bias_f: m_bias_fs) {
       delete bias_f;
  }
}

void SystemBiases::setup_distance_bias(InputParameters& params) {

    vector<OrderParam*> order_params {m_system_order_params.get_distance_params()};
    for (auto order_param: order_params) {

        // For now only use the linear step function
        BiasFunction* bias_f;
        bias_f = new LinearStepBiasFunction {*order_param, params.m_min_dist,
                params.m_max_dist, params.m_max_bias};
        m_bias_fs.push_back(bias_f);

        vector<Domain*> domains {order_param->get_depending_domains()};
        for (auto domain: domains) {
            add_bias_f_dependency(*domain, bias_f);
        }
    }
}

vector<BiasFunction*> SystemBiases::get_dependent_bias_fs(Domain& domain) {
    pair<int, int> key {domain.m_c, domain.m_d};
    return m_domain_to_bias_fs[key];
}

void SystemBiases::add_bias_f_dependency(Domain& domain, BiasFunction* bias_f) {
    int c_i {domain.m_c};
    int d_i {domain.m_d};
    pair<int, int> key {c_i, d_i};
    m_domain_to_bias_fs[key].push_back(bias_f);
}

//void SystemBiases::remove_domain_dependencies(Domain& domain) {
//}

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

double SystemBiases::calc_bias() {
    double total_bias {0};
    for (auto bias_f: m_bias_fs) {
        double bias {bias_f->update_bias()};
        total_bias += bias;
    }

    return total_bias * m_bias_mult;
}

void SystemBiases::update_bias_mult(double bias_mult) {
    m_bias_mult = bias_mult;
}

double SystemBiases::calc_one_domain(Domain& domain) {
    double bias_diff {0};

    // Get dependent bias functions and update
    m_system_order_params.update_one_domain(domain);
    vector<BiasFunction*> bias_fs {get_dependent_bias_fs(domain)};
    for (auto bias_f: bias_fs) {
        double prev_bias {bias_f->get_bias()};
        double bias {bias_f->update_bias()};
        bias_diff += bias - prev_bias;
    }

    return bias_diff;
}

double SystemBiases::check_one_domain(Domain& domain, VectorThree new_pos,
        VectorThree new_ore, Occupancy new_state) {
    double bias_diff {0};

    // Get dependent bias functions
    vector<BiasFunction*> bias_fs {get_dependent_bias_fs(domain)};

    // Check bias bias values without internal update
    for (auto bias_f: bias_fs) {
        double prev_bias {bias_f->get_bias()};
        double bias {bias_f->check_bias(domain, new_pos, new_ore, new_state)};
        bias_diff += bias - prev_bias;
    }

    return bias_diff;
}
