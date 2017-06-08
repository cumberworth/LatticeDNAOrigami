// order_params.cpp

#include<numeric>

#include "utility.h"
#include "order_params.h"

using namespace Utility;
using namespace OrderParams;

using std::abs;

string OrderParam::get_label() {
    return m_label;
}

DistOrderParam::DistOrderParam(Domain& domain_1, Domain& domain_2) :
        m_domain_1 {domain_1},
        m_domain_2 {domain_2} {
    calc_param();
}

int DistOrderParam::calc_param() {
    if (m_domain_1.m_state != Occupancy::unassigned and
            m_domain_2.m_state != Occupancy::unassigned) {

        m_defined = true;
        VectorThree diff_vec {m_domain_2.m_pos - m_domain_1.m_pos};
        m_param = diff_vec.abssum();
        m_checked_param = m_param;
    }
    else {
        m_defined = false;
    }

    return m_param;
}

int DistOrderParam::check_param(Domain& domain, VectorThree new_pos, VectorThree,
        Occupancy state) {
    Domain* unmodded_domain {&m_domain_1};
    if (domain.m_d == m_domain_1.m_d) {
        unmodded_domain = &m_domain_2;
    }
    if (unmodded_domain->m_state != Occupancy::unassigned and
            state != Occupancy::unassigned) {
        m_checked_param = (new_pos - unmodded_domain->m_pos).abssum();
        // Will be defined for configurations where that domain is unassigned
        // Consider using a check defined variable instead
        m_defined = true;
    }
    else {
        m_defined = false;
    }

    return m_checked_param;
}

int DistOrderParam::get_param() {
    return m_param;
}

int DistOrderParam::get_checked_param() {
    return m_checked_param;
}

string DistOrderParam::get_label() {
    return m_label;
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

DistSumOrderParam::DistSumOrderParam(vector<DistOrderParam*> dist_params) :
        m_dist_params {dist_params} {
    calc_param();
}

int DistSumOrderParam::calc_param() {
    int dist_sum {0};
    m_defined = true;
    for (auto param: m_dist_params) {

        // Will be defined if one of the domains is unassigned but being checked
        // Obvoiusly this comment is prone to being out-of-date
        if (param->defined()) {
            dist_sum += param->get_param();
        }
        else {
            m_defined = false;
            dist_sum = m_param;
            break;
        }
    }
    m_param = dist_sum;
    m_checked_param = dist_sum;

    return m_param;
}

int DistSumOrderParam::check_param(Domain&, VectorThree, VectorThree,
        Occupancy) {
    int dist_sum {0};
    m_defined = true;
    for (auto param: m_dist_params) {

        // Will be defined if one of the domains is unassigned but being checked
        // Obvoiusly this comment is prone to being out-of-date
        if (param->defined()) {
            dist_sum += param->get_checked_param();
        }
        else {
            m_defined = false;
            dist_sum = m_param;
            break;
        }
    }
    m_checked_param = dist_sum;

    return dist_sum;
}

int DistSumOrderParam::get_param() {
    return m_param;
}

int DistSumOrderParam::get_checked_param() {
    return m_checked_param;
}

string DistSumOrderParam::get_label() {
    return m_label;
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
    return m_defined;
}

NumStaplesOrderParam::NumStaplesOrderParam(OrigamiSystem& origami) :
        m_origami {origami} {
    calc_param();
}

int NumStaplesOrderParam::calc_param() {

    // This is only defined for fully set configurations
    //HACK
    /*if (m_origami.configuration_fully_set()) {
        m_defined = true;
        m_param = m_origami.num_staples();
        m_checked_param = m_param;
    }
    else {
        m_defined = false;
    }*/
    m_param = m_origami.num_staples();

    return m_param;
}

int NumStaplesOrderParam::check_param(Domain&, VectorThree , VectorThree,
        Occupancy occ) {

    // This is only defined for fully set configurations
    if (m_origami.num_unassigned_domains() == 1 and
            occ != Occupancy::unassigned) {
        m_defined = true;
        m_checked_param = m_origami.num_staples();
    }
    else {
        m_defined = false;
    }
    return m_checked_param;
}

int NumStaplesOrderParam::check_delete_chain(int) {
    // Actually just want it to be defined now, only ever delete one chain at
    // a time and never unassign any other domains
    m_defined = true;
    m_checked_param = m_origami.num_staples() - 1;

    return m_checked_param;
}

int NumStaplesOrderParam::get_param() {
    return m_param;
}

int NumStaplesOrderParam::get_checked_param() {
    return m_checked_param;
}

string NumStaplesOrderParam::get_label() {
    return m_label;
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
    return m_defined;
}

NumBoundDomainPairsOrderParam::NumBoundDomainPairsOrderParam(
        OrigamiSystem& origami) :
        m_origami {origami} {
    calc_param();
}

int NumBoundDomainPairsOrderParam::calc_param() {

    // This is only defined for fully set configurations
    //HACK
    /*if (m_origami.configuration_fully_set()) {
        m_param = m_origami.num_fully_bound_domain_pairs();
        m_checked_param = m_param;
        m_defined = true;
    }
    else {
        m_defined = false;
    }*/
    m_param = m_origami.num_fully_bound_domain_pairs();

    return m_param;
}

int NumBoundDomainPairsOrderParam::check_param(Domain&, VectorThree , VectorThree,
        Occupancy state) {
    
    // This is only defined for fully set configurations
    if (m_origami.num_unassigned_domains() == 1) {
        if (state != Occupancy::unassigned) {
            m_defined = true;
            int num_domain_pairs = m_origami.num_fully_bound_domain_pairs();
            if (state == Occupancy::bound) {
                m_checked_param = num_domain_pairs + 1;
            }
            else {
                m_checked_param = num_domain_pairs;
            }
        }
        else {
            m_defined = false;
        }
    }

    return m_checked_param;
}

int NumBoundDomainPairsOrderParam::check_delete_chain(int) {
    // Actually just want it to be defined now, only ever delete one chain at
    // a time and never unassign any other domains
    m_defined = true;
    m_checked_param = m_origami.num_fully_bound_domain_pairs();

    return m_checked_param;
}

int NumBoundDomainPairsOrderParam::get_param() {
    return m_param;
}

int NumBoundDomainPairsOrderParam::get_checked_param() {
    return m_checked_param;
}

string NumBoundDomainPairsOrderParam::get_label() {
    return m_label;
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
    return m_defined;
}

SystemOrderParams::SystemOrderParams(InputParameters& params,
        OrigamiSystem& origami) :
        m_origami {origami},
        m_num_staples {origami},
        m_num_bound_domains {origami} {

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
        m_dist_sums.push_back(dist_sum);
        m_order_params.push_back(dist_sum);
    }

    m_order_params.push_back(&m_num_staples);
    m_order_params.push_back(&m_num_bound_domains);
}

SystemOrderParams::~SystemOrderParams() {
    for (auto op: m_dists) {
        delete op;
    }
    for (auto op: m_dist_sums) {
        delete op;
    }
}

vector<DistOrderParam*> SystemOrderParams::get_distance_params() {
    // For now only distance order params, but will need to change eventually
    return m_dists;
}

vector<DistSumOrderParam*> SystemOrderParams::get_dist_sums() {
    return m_dist_sums;
}

NumBoundDomainPairsOrderParam& SystemOrderParams::get_num_bound_domains() {
    return m_num_bound_domains;
}

NumStaplesOrderParam& SystemOrderParams::get_num_staples() {
    return m_num_staples;
}

vector<OrderParam*> SystemOrderParams::get_order_params() {
    return m_order_params;
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
        DistOrderParam* order_param;
        order_param = new DistOrderParam {domain_1, domain_2};
        m_dists.push_back(order_param);
        m_order_params.push_back(order_param);
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

vector<OrderParam*> SystemOrderParams::get_dependent_dists(
        Domain& domain) {
    pair<int, int> key {domain.m_c, domain.m_d};
    return m_domain_to_order_params[key];
}

void SystemOrderParams::update_params() {
    // HACK
    for (auto order_param: m_order_params) {
        order_param->calc_param();
    }
}

void SystemOrderParams::update_params_deletion_case() {
    // SUPER HACK
    for (auto order_param: m_order_params) {
        order_param->calc_param();
    }
    // SUPER SUPER HACK
    m_num_staples.set_param(m_num_staples.get_param() - 1);
}

void SystemOrderParams::update_one_domain(Domain& domain) {
    // Do distance ones, then sum, and also the num bound domains and staples.
    // There will a seperate vector for each order parameter type; I just iterate
    // through each. I can always make it more general in the future so that I don't
    // have to manually put the order parameters a certain way to ensure that
    // ones that depend on others have had those ones updated first.

    // Get distances parameters that are dependent on the domain
    vector<OrderParam*> dependent_dists {get_dependent_dists(domain)};
    for (auto dist: dependent_dists) {
        dist->calc_param();
    }

    // Update order params dependent on dists
    for (auto dist_sum: m_dist_sums) {
        if (dist_sum->dependent_on(domain)) {
            dist_sum->calc_param();
        }
    }

    // For now always calculate as it's just getting a pre-calculated value
    m_num_staples.calc_param();
    m_num_bound_domains.calc_param();
}

void SystemOrderParams::update_delete_chain(int) {

    // For now always only these depend on deletion of staple chain
    m_num_staples.calc_param();
    m_num_bound_domains.calc_param();
}

void SystemOrderParams::check_one_domain(Domain& domain, VectorThree pos,
        VectorThree ore, Occupancy state) {
    // Get distance parameters that are dependent on the domain
    vector<OrderParam*> dependent_dists {get_dependent_dists(domain)};
    for (auto dist: dependent_dists) {
        dist->check_param(domain, pos, ore, state);
    }

    // Update order params dependent on dists
    for (auto dist_sum: m_dist_sums) {
        if (dist_sum->dependent_on(domain)) {
            dist_sum->check_param(domain, pos, ore, state);
        }
    }

    // For now always calculate as it's just getting a pre-calculated value
    m_num_staples.check_param(domain, pos, ore, state);
    m_num_bound_domains.check_param(domain, pos, ore, state);
}

void SystemOrderParams::check_delete_chain(int c_i) {

    // For now always only these depend on deletion of staple chain
    m_num_staples.check_delete_chain(c_i);
    m_num_bound_domains.check_delete_chain(c_i);
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
    //HACK
    //if (m_order_param.defined()) {
        new_bias = calc_bias(param);
    //}
    //else {
    //    new_bias = 0;
    //}
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

bool SquareWellBiasFunction::dependent_on(OrderParam& order_param) {
    bool dependent {false};
    if (&order_param == &m_order_param) {
        dependent = true;
    }

    return dependent;
}

double SquareWellBiasFunction::get_bias() {
    return m_bias;
}

GridBiasFunction::GridBiasFunction() {
}

void GridBiasFunction::set_order_params(vector<OrderParam*> order_params) {
    m_order_params = order_params;
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
    for (auto order_param: m_order_params) {

        //HACK
        // No bias if order parameters not all defined
        //if (not order_param->defined()) {
        //    m_bias = 0;
        //    return 0;
        //}
        key.push_back(order_param->get_param());
    }
    double bias {calc_bias(key)};
    m_bias = bias;

    return bias;
}

double GridBiasFunction::check_bias() {
    vector<int> key {};
    for (auto order_param: m_order_params) {

        // No bias if order parameters not all defined
        if (not order_param->defined()) {
            return 0;
        }
        key.push_back(order_param->get_checked_param());
    }
    double checked_bias {calc_bias(key)};

    return checked_bias;
}

bool GridBiasFunction::dependent_on(OrderParam& order_param) {
    bool dependent {false};
    for (auto dependent_param: m_order_params) {
        if (&order_param == dependent_param) {
            dependent = true;
        }
        else {
            dependent = false;
        }
    }

    return dependent;
}

double GridBiasFunction::get_bias() {
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
    OrderParam* op1 {&m_system_order_params.get_num_staples()};
    OrderParam* op2 {&m_system_order_params.get_num_bound_domains()};
    if (params.m_square_well_bias) {
       add_square_well_bias(op1, params.m_min_well_param, params.m_max_well_param,
             params.m_well_bias, params.m_outside_bias);
       add_square_well_bias(op2, params.m_min_well_param, params.m_max_well_param,
             params.m_well_bias, params.m_outside_bias);
       //add_square_well_bias(op, 2, 4, 0, 99);
    //op = &m_system_order_params.get_num_staples();
       //add_square_well_bias(op, 2, 3, 0, 99);
    }
    m_bias_fs.push_back(&m_grid_bias_f);

    // Update total bias
    for (auto bias_f: m_bias_fs) {
        m_total_bias += bias_f->get_bias();
    }
}

SystemBiases::~SystemBiases() {
    for (auto bias_f: m_dist_bias_fs) {
        delete bias_f;
    }
    for (auto bias_f: m_well_bias_fs) {
        delete bias_f;
    }
}

void SystemBiases::setup_distance_bias(InputParameters& params) {

    vector<DistOrderParam*> order_params {
            m_system_order_params.get_distance_params()};
    for (auto order_param: order_params) {

        // For now only use the linear step function
        LinearStepBiasFunction* bias_f;
        bias_f = new LinearStepBiasFunction {*order_param, params.m_min_dist,
                params.m_max_dist, params.m_max_bias};
        m_dist_bias_fs.push_back(bias_f);
        m_bias_fs.push_back(bias_f);

        vector<Domain*> domains {order_param->get_depending_domains()};
        for (auto domain: domains) {
            add_bias_f_dependency(*domain, bias_f);
        }
    }
}

vector<BiasFunction*> SystemBiases::get_dependent_dist_restraints(Domain& domain) {
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
    m_system_order_params.update_params();
    double total_bias {0};
    for (auto bias_f: m_bias_fs) {
        double bias {bias_f->update_bias()};
        total_bias += bias;
    }
    m_total_bias = total_bias;

    return total_bias * m_bias_mult;
}

double SystemBiases::calc_bias_delete_case() {
    m_system_order_params.update_params_deletion_case();
    double total_bias {0};
    for (auto bias_f: m_bias_fs) {
        double bias {bias_f->update_bias()};
        total_bias += bias;
    }
    m_total_bias = total_bias;

    return total_bias * m_bias_mult;
}

double SystemBiases::get_bias() {
    return m_total_bias * m_bias_mult;
}

void SystemBiases::update_bias_mult(double bias_mult) {
    m_bias_mult = bias_mult;
}

double SystemBiases::calc_one_domain(Domain& domain) {
    double bias_diff {0};

    // Get dependent bias functions and update
    m_system_order_params.update_one_domain(domain);
    vector<BiasFunction*> dist_bias_fs {get_dependent_dist_restraints(domain)};
    for (auto dist_bias_f: dist_bias_fs) {
        double prev_bias {dist_bias_f->get_bias()};
        double new_bias {dist_bias_f->update_bias()};
        bias_diff += new_bias - prev_bias;
    }

    // Square well bias fs
    for (auto well_bias_f: m_well_bias_fs) {
        double prev_bias {well_bias_f->get_bias()};
        double new_bias {well_bias_f->update_bias()};
        bias_diff += new_bias - prev_bias;
    }

    // Grid bias f
    double prev_bias {m_grid_bias_f.get_bias()};
    double new_bias {m_grid_bias_f.update_bias()};
    bias_diff += new_bias - prev_bias;
    m_total_bias += bias_diff;

    return bias_diff * m_bias_mult;
}

double SystemBiases::calc_delete_chain(int c_i) {
    m_system_order_params.update_delete_chain(c_i);

    // For now assuming no staple distance constraints
    double prev_bias {m_grid_bias_f.get_bias()};
    double new_bias {m_grid_bias_f.update_bias()};
    double bias_diff {new_bias - prev_bias};


    // Square well
    for (auto well_bias_f: m_well_bias_fs) {
        double prev_bias {well_bias_f->get_bias()};
        double new_bias {well_bias_f->update_bias()};
        bias_diff += new_bias - prev_bias;
    }

    m_total_bias += bias_diff;

    return bias_diff * m_bias_mult;
}

double SystemBiases::check_one_domain(Domain& domain, VectorThree pos,
        VectorThree ore, Occupancy state) {
    double bias_diff {0};
    m_system_order_params.check_one_domain(domain, pos, ore, state);

    // Get dependent bias functions
    vector<BiasFunction*> dist_bias_fs {get_dependent_dist_restraints(domain)};

    // Check bias bias values without internal update
    for (auto dist_bias_f: dist_bias_fs) {
        double prev_bias {dist_bias_f->get_bias()};
        double new_bias {dist_bias_f->check_bias()};
        bias_diff += new_bias - prev_bias;
    }

    // Square well
    for (auto well_bias_f: m_well_bias_fs) {
        double prev_bias {well_bias_f->get_bias()};
        double new_bias {well_bias_f->check_bias()};
        bias_diff += new_bias - prev_bias;
    }

    // Grid
    double prev_bias {m_grid_bias_f.get_bias()};
    double new_bias {m_grid_bias_f.check_bias()};
    bias_diff += new_bias - prev_bias;

    return bias_diff * m_bias_mult;
}

double SystemBiases::check_delete_chain(int c_i) {
    m_system_order_params.check_delete_chain(c_i);

    // For now only effects these biases
    double prev_bias {m_grid_bias_f.get_bias()};
    double new_bias {m_grid_bias_f.check_bias()};
    double bias_diff {new_bias - prev_bias};

    // Square well
    for (auto well_bias_f: m_well_bias_fs) {
        double prev_bias {well_bias_f->get_bias()};
        double new_bias {well_bias_f->check_bias()};
        bias_diff += new_bias - prev_bias;
    }

    return bias_diff * m_bias_mult;
}

GridBiasFunction* SystemBiases::get_grid_bias() {
    return &m_grid_bias_f;
}

void SystemBiases::add_square_well_bias(OrderParam* op, int well_min, int well_max,
        double well_bias, double outside_bias) {
    SquareWellBiasFunction* well_bias_f;
    well_bias_f = new SquareWellBiasFunction {*op,
        well_min, well_max, well_bias, outside_bias};
    m_well_bias_fs.push_back(well_bias_f);
    m_bias_fs.push_back(well_bias_f);

    // Update total bias
    m_total_bias += well_bias_f->get_bias();
}
