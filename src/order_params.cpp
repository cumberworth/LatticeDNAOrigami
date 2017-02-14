// order_params.cpp

#include "utility.h"
#include "order_params.h"

using namespace Utility;
using namespace OrderParams;

DistOrderParam::DistOrderParam(Domain& domain_1, Domain& domain_2) :
        m_domain_1 {domain_1},
        m_domain_2 {domain_2} {
}

int DistOrderParam::calc_param() {
    VectorThree diff_vec {m_domain_2.m_pos - m_domain_1.m_pos};
    return diff_vec.sum();
}

LinearStepBiasFunction::LinearStepBiasFunction(OrderParam& order_param,
        int min_param, int max_param, double max_bias) :
        m_order_param {order_param},
        m_min_param {min_param},
        m_max_param {max_param},
        m_max_bias {max_bias} {

    m_slope = m_max_bias / (max_param - min_param);
}

LinearStepBiasFunction::~LinearStepBiasFunction() {
    delete &m_order_param;
}

double LinearStepBiasFunction::calc_bias() {
    double bias {0};
    int param {m_order_param.calc_param()};
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

SystemBias::SystemBias(InputParameters& params, OrigamiSystem& origami) {
    // For now everything is assumed to be distance order parameters with linear
    // function biases with a single min and max param
    m_bias_mult = params.m_bias_mult;

    // Extract domain ids and store domain pairs
    vector<pair<Domain&, Domain&>> domain_pairs {};
    size_t num_restrained_domains {params.m_restraint_pairs.size()};
    for (size_t pair_i {0}; pair_i != num_restrained_domains; pair_i += 2) {
        int domain_1_i {params.m_restraint_pairs[pair_i]};
        Domain& domain_1 {*origami.get_domain(0, domain_1_i)}; 
        int domain_2_i {params.m_restraint_pairs[pair_i + 1]};
        Domain& domain_2 {*origami.get_domain(0, domain_2_i)};
        domain_pairs.push_back({domain_1, domain_2});
    }

    // Setup up order parameters and biases
    for (size_t pair_i {0}; pair_i != domain_pairs.size(); pair_i++) {
        Domain& domain_1 {domain_pairs[pair_i].first};
        Domain& domain_2 {domain_pairs[pair_i].second};
        OrderParam* op;
        op = new DistOrderParam {domain_1, domain_2};
        BiasFunction* bias_f;
        bias_f = new LinearStepBiasFunction {*op, params.m_min_dist, params.m_max_dist,
            params.m_max_bias};
        m_biases.push_back(bias_f);
    }
}

SystemBias::~SystemBias() {
    for (auto bias_f: m_biases) {
        delete bias_f;
    }
}

double SystemBias::calc_bias() {
    double bias {0};
    for (auto bias_f: m_biases) {
        bias += bias_f->calc_bias();
    }

    return bias * m_bias_mult;
}

void SystemBias::update_bias_mult(double bias_mult) {
    m_bias_mult = bias_mult;
}
