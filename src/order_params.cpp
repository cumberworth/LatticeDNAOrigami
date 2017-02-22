// order_params.cpp

#include "utility.h"
#include "order_params.h"

using namespace Utility;
using namespace OrderParams;

DistOrderParam::DistOrderParam(Domain& domain_1, Domain& domain_2) :
        m_domain_1 {domain_1},
        m_domain_2 {domain_2} {
    calc_param();
}

int DistOrderParam::calc_param() {
    int dist;
        VectorThree diff_vec {m_domain_2.m_pos - m_domain_1.m_pos};
        dist = diff_vec.sum();
    if (m_domain_1.m_state != Occupancy::unassigned and
            m_domain_2.m_state != Occupancy::unbound) {
        m_defined = true;
    }
    else {
        m_defined = false;
    }

    return dist;
}

int DistOrderParam::calc_param(Domain& domain, VectorThree new_pos, VectorThree,
        Occupancy) {
    VectorThree pos_2;
    if (domain.m_d == m_domain_1.m_d) {
        pos_2 = m_domain_2.m_pos;
    }
    else {
        pos_2 = m_domain_1.m_pos;
    }

    return (new_pos - pos_2).sum();
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

bool DistOrderParam::defined() {
    return m_defined;
}

LinearStepBiasFunction::LinearStepBiasFunction(OrderParam& order_param,
        int min_param, int max_param, double max_bias) :
        m_order_param {order_param},
        m_min_param {min_param},
        m_max_param {max_param},
        m_max_bias {max_bias} {

    m_slope = m_max_bias / (max_param - min_param);
    calc_bias();
}

LinearStepBiasFunction::~LinearStepBiasFunction() {
    delete &m_order_param;
}

double LinearStepBiasFunction::calc_bias() {
    double bias;
    if (m_order_param.defined()) {
        int param {m_order_param.get_param()};
        if (param <= m_min_param) {
            bias = 0;
        }
        else if (param > m_min_param and param <= m_max_param) {
            bias = m_slope * (param - m_min_param);
        }
        else {
            bias = m_max_bias;
        }
    }
    else {
        bias = 0;
    }
    
    return bias;
}

double LinearStepBiasFunction::update_bias() {
    double bias {calc_bias()};
    m_bias = bias;

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

SystemBias::SystemBias(InputParameters& params, OrigamiSystem& origami) :
        m_origami {origami},
        m_bias_mult {params.m_bias_mult} {

    // Setup dependency tables
    for (auto chain: origami.get_chains()) {
        for (auto domain: chain) {
            pair<int, int> key {domain->m_c, domain->m_d};
            m_domain_to_order_params[key] = {};
            m_domain_to_bias_fs[key] = {};
        }
    }

    // Setup each type of order parameter and bias function
    if (params.m_distance_bias) {
        setup_distance_bias(params);
    }
}

SystemBias::~SystemBias() {
    for (auto bias_f: m_bias_fs) {
       delete bias_f;
  }
}

void SystemBias::setup_distance_bias(InputParameters& params) {

    // Extract domain ids and store domain pairs
    // Consider this a part of the parser module
    vector<pair<Domain&, Domain&>> domain_pairs {};
    size_t num_restrained_domains {params.m_restraint_pairs.size()};
    for (size_t pair_i {0}; pair_i != num_restrained_domains; pair_i += 2) {
        int domain_1_i {params.m_restraint_pairs[pair_i]};
        Domain& domain_1 {*m_origami.get_domain(0, domain_1_i)}; 
        int domain_2_i {params.m_restraint_pairs[pair_i + 1]};
        Domain& domain_2 {*m_origami.get_domain(0, domain_2_i)};
        domain_pairs.push_back({domain_1, domain_2});
    }

    // Setup and store order parameters and biases
    for (size_t pair_i {0}; pair_i != domain_pairs.size(); pair_i++) {
        Domain& domain_1 {domain_pairs[pair_i].first};
        Domain& domain_2 {domain_pairs[pair_i].second};

        OrderParam* order_param;
        order_param = new DistOrderParam {domain_1, domain_2};
        m_order_params.push_back(order_param);
        add_order_param_dependency(domain_1, order_param);

        // For now only use the linear step function
        BiasFunction* bias_f;
        bias_f = new LinearStepBiasFunction {*order_param, params.m_min_dist,
                params.m_max_dist, params.m_max_bias};
        m_bias_fs.push_back(bias_f);
        add_bias_f_dependency(domain_1, bias_f);
    }
}

vector<OrderParam*> SystemBias::get_dependent_order_params(Domain& domain) {
    pair<int, int> key {domain.m_c, domain.m_d};
    return m_domain_to_order_params[key];
}

vector<BiasFunction*> SystemBias::get_dependent_bias_fs(Domain& domain) {
    pair<int, int> key {domain.m_c, domain.m_d};
    return m_domain_to_bias_fs[key];
}

void SystemBias::add_order_param_dependency(Domain& domain, OrderParam* order_param) {
    int c_i {domain.m_c};
    int d_i {domain.m_d};
    pair<int, int> key {c_i, d_i};
    m_domain_to_order_params[key].push_back(order_param);
}

void SystemBias::add_bias_f_dependency(Domain& domain, BiasFunction* bias_f) {
    int c_i {domain.m_c};
    int d_i {domain.m_d};
    pair<int, int> key {c_i, d_i};
    m_domain_to_bias_fs[key].push_back(bias_f);
}

//void SystemBias::remove_domain_dependencies(Domain& domain) {
//}

//void SystemBias::add_chain(vector<Domain*> chain) {
    // For each domain
    // If biases that include all chains of type x
        // Create new order parameters and bias functions
        // Add dependencies
//}

//void SystemBias::remove_chain(vector<Domain*> chain) {
    // For each domain
    // Find order params and bias fs dependent on domain
    // Remove dependecies
    // Iterate through order params and bias fs
    //     delete those with 0 dependencies
//}

double SystemBias::calc_bias() {
    double total_bias {0};
    for (auto bias_f: m_bias_fs) {
        double bias {bias_f->calc_bias()};
        total_bias += bias;
    }

    return total_bias * m_bias_mult;
}

void SystemBias::update_bias_mult(double bias_mult) {
    m_bias_mult = bias_mult;
}

double SystemBias::calc_one_domain(Domain& domain) {
    double bias_diff {0};

    // Get dependent order params and bias functions and update
    vector<OrderParam*> order_params {get_dependent_order_params(domain)};
    vector<BiasFunction*> bias_fs {get_dependent_bias_fs(domain)};
    for (auto order_param: order_params) {
        order_param->calc_param();
    }
    for (auto bias_f: bias_fs) {
        double prev_bias {bias_f->get_bias()};
        double bias {bias_f->update_bias()};
        bias_diff += bias - prev_bias;
    }

    return bias_diff;
}

double SystemBias::check_one_domain(Domain& domain, VectorThree new_pos,
        VectorThree new_ore, Occupancy new_state) {
    double bias_diff {0};

    // Get dependent order params and bias functions
    vector<OrderParam*> order_params {get_dependent_order_params(domain)};
    vector<BiasFunction*> bias_fs {get_dependent_bias_fs(domain)};

    // Calculate order parameters with trial values
    for (auto order_param: order_params) {
        order_param->calc_param(domain, new_pos, new_ore, new_state);
    }

    // Calculate bias function without updating current bias
    for (auto bias_f: bias_fs) {
        double prev_bias {bias_f->get_bias()};
        double bias {bias_f->calc_bias()};
        bias_diff += bias - prev_bias;
    }

    // Revert order parameters to original state. This is probably not
    // necessary; I will always finish with one domain before moving
    // on to another, so it would always end up being set correctly
    for (auto order_param: order_params) {
        order_param->calc_param();
    }

    return bias_diff;
}
