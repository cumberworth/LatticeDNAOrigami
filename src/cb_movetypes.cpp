// cb_movetypes.cpp

#include "utility.h"
#include "movetypes.h"
#include "cb_movetypes.h"

using namespace Movetypes;
using namespace CBMovetypes;
using namespace Utility;

void CBMCMovetype::reset_internal() {
    MCMovetype::reset_internal();
    m_bias = 1;
    m_new_bias = 1;
    m_new_modifier = 1;
    m_regrow_old = false;
    m_old_pos.clear();
    m_old_ore.clear();
}

void CBMCMovetype::calc_biases(
        Domain& domain,
        VectorThree p_prev,

        // Pairs of position and orientations (if bound, orientation
        // vector, otherwise the number of possible orientations (6)
        vector<pair<VectorThree, VectorThree>>& configs,

        // List of associated boltzmann factors
        vector<double>& bfactors) {

    // Iterate through all possible new positions
    for (auto v: vectors) {

        // Trial position vector
        VectorThree p_new {p_prev + v};

        // Check energies of each configuration
        switch (m_origami_system.position_occupancy(p_new)) {
            case Occupancy::bound:
            case Occupancy::misbound:
                continue;
            case Occupancy::unbound: {
                Domain* unbound_domain {m_origami_system.unbound_domain_at(p_new)};
                //bool domains_complementary {m_origami_system.check_domains_complementary(
                //        domain, *unbound_domain)};
                VectorThree o_new;
                //double bfactor;
                //if (domains_complementary) {
                //    o_new = -unbound_domain->m_ore;
                //    bfactor = 1;
                //}
                //else {
                //    o_new = {0, 0, 0};
                //    bfactor = 6;
                //}
                double bfactor {1};
                o_new = -unbound_domain->m_ore;
                double delta_e {m_origami_system.check_domain_constraints(
                        domain, p_new, o_new)};
                if (not m_origami_system.m_constraints_violated) {
                    configs.push_back({p_new, o_new});
                    bfactor *= exp(-delta_e);
                    bfactors.push_back(bfactor);
                }
                else {
                    m_origami_system.m_constraints_violated = false;
                }
                break;
            }
            case Occupancy::unassigned:
                configs.push_back({p_new, {0, 0, 0}});

                // Biases can be position dependent, but will be orientation
                // independent
                double delta_e {m_origami_system.check_domain_constraints(
                        domain, p_new, {1, 0, 0})};
                bfactors.push_back(6 * exp(-delta_e));
        }
    }
}

void CBMCMovetype::select_and_set_config(vector<Domain*> domains, int i) {
    Domain* domain {domains[i]};
    Domain* prev_domain {domains[i - 1]};
    VectorThree p_prev {prev_domain->m_pos};

    // Calculate weights
    vector<pair<VectorThree, VectorThree>> configs {};
    vector<double> bfactors {};
    calc_biases(*domain, p_prev, configs, bfactors);
    vector<double> weights {calc_bias(bfactors, domain, configs, p_prev, domains)};
    if (m_rejected) {
        return;
    }

    if (not m_regrow_old) {
        select_and_set_new_config(*domain, weights, configs);
    }
    else {
        select_and_set_old_config(*domain);
    }

    // Reversion list update
    pair<int, int> key {domain->m_c, domain->m_d};
    m_assigned_domains.push_back(key);
}

void CBMCMovetype::select_and_set_new_config(
        Domain& domain,
        vector<double> weights,
        vector<pair<VectorThree, VectorThree>> configs) {

    // Select config based on weights
    pair<VectorThree, VectorThree> config;
    double cum_prob {0};
    double random_real {m_random_gens.uniform_real()};
    for (size_t i {0}; i != weights.size(); i++) {
        cum_prob += weights[i];
        if (random_real < cum_prob) {
            config = configs[i];
            break;
        }
    }

    // Set config on domain
    VectorThree p_new {config.first};
    VectorThree o_new {config.second};
    if (o_new == VectorThree {0, 0, 0}) {
        o_new = select_random_orientation();
        m_origami_system.set_checked_domain_config(domain, p_new, o_new);
    }
    else {
        m_origami_system.set_checked_domain_config(domain, p_new, o_new);
    }
}

void CBMCMovetype::select_and_set_old_config(Domain& domain) {
    pair<int, int> key {domain.m_c, domain.m_d};
    VectorThree p_old {m_old_pos[key]};
    VectorThree o_old {m_old_ore[key]};
    m_origami_system.set_checked_domain_config(domain, p_old, o_old);
}

//double CBMCMovetype::set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old) {
//    bool domains_complementary {m_origami_system.check_domains_complementary(
//            growth_domain_new, growth_domain_old)};
//    VectorThree o_new;
//    double bias {1};
//    if (domains_complementary) {
//        o_new = -growth_domain_old.m_ore;
//    }
//    else {
//        o_new = select_random_orientation();
//        bias *= 6;
//    }
//    double delta_e {0};
//    delta_e += m_origami_system.set_domain_config(growth_domain_new,
//            growth_domain_old.m_pos, o_new);
//    if (m_origami_system.m_constraints_violated) {
//        m_rejected = true;
//    }
//    else {
//        bias *= exp(-delta_e);
//        m_bias *= bias;
//        pair<int, int> key {growth_domain_new.m_c, growth_domain_new.m_d};
//        m_assigned_domains.push_back(key);
//    }
//
//    return delta_e;
//}

double CBMCMovetype::set_old_growth_point(Domain& growth_domain_new, Domain& growth_domain_old) {
    pair<int, int> key {growth_domain_new.m_c, growth_domain_new.m_d};
    VectorThree o_old {m_old_ore[key]};
    double delta_e {0};
    delta_e += m_origami_system.set_checked_domain_config(growth_domain_new,
            growth_domain_old.m_pos, o_old);
    //bool domains_complementary {m_origami_system.check_domains_complementary(
    //        growth_domain_new, growth_domain_old)};
    //double bias {1};
    //if (not domains_complementary) {
    //    bias *= 6;
    //}
    //bias *= exp(-delta_e);
    m_bias *= exp(-delta_e);
    m_assigned_domains.push_back(key);

    return delta_e;
}

bool CBMCMovetype::test_cb_acceptance() {
    long double ratio {m_new_bias / m_bias};
    bool accepted;
    if (test_acceptance(ratio)) {
        reset_origami();
        accepted = true;
    }
    else {
        m_modified_domains.clear();
        m_assigned_domains.clear();
        accepted = false;
    }
    return accepted;
}

void CBMCMovetype::unassign_domains(vector<Domain*> domains) {
    for (auto domain: domains) {
        pair<int, int> key {domain->m_c, domain->m_d};
        m_prev_pos[key] = domain->m_pos;
        m_prev_ore[key] = domain->m_ore;
        m_modified_domains.push_back(key);

        // Would be double counting to include delta_e in bias
        m_origami_system.unassign_domain(*domain);
    }
}

void CBMCMovetype::setup_for_regrow_old() {
    // Save relevant state variables for acceptance testing and resetting
    m_regrow_old = true;
    m_new_bias = m_bias;
    m_new_modifier = m_modifier;
    m_bias = 1;
    m_modified_domains.clear();
    m_assigned_domains.clear();
    m_old_pos = m_prev_pos;
    m_old_ore = m_prev_ore;
}

vector<pair<Domain*, Domain*>> CBMCMovetype::find_bound_domains(
        vector<Domain*> selected_chain) {

    vector<pair<Domain*, Domain*>> bound_domains {};
    for (auto domain: selected_chain) {
        // shouldn't this be only non-self binding (only would effect staple size > 2)
        if (domain->m_bound_domain != nullptr) {

            // New domain, old domain
            bound_domains.push_back({domain, domain->m_bound_domain});
        }
    }

    if (bound_domains.empty()) {
        cout << "d\n";
        throw OrigamiMisuse {};
    }

    return bound_domains;
}

pair<Domain*, Domain*> CBMCMovetype::select_old_growthpoint(
        vector<pair<Domain*, Domain*>> bound_domains) {

    int bound_domain_index {m_random_gens.uniform_int(0, bound_domains.size() - 1)};
    Domain* growth_domain_new {bound_domains[bound_domain_index].first};
    Domain* growth_domain_old {bound_domains[bound_domain_index].second};
    return {growth_domain_new, growth_domain_old};
}

void CBMCMovetype::update_bias(int sign) {
    double total_bias {m_system_bias.calc_bias()};
    m_bias *= exp(-sign*total_bias);
}

bool CBStapleRegrowthMCMovetype::attempt_move() {
    bool accepted;

    // No staples to regrow
    if (m_origami_system.num_staples() == 0) {
        accepted = false;
        return accepted;
    }

    // Select a staple to regrow
    int c_i_index {m_random_gens.uniform_int(1, m_origami_system.num_staples())};
    vector<Domain*> selected_chain {m_origami_system.get_chains()[c_i_index]};

    // Reject if staple is connector
    if (staple_is_connector(selected_chain)) {
        accepted = false;
        return accepted;
    }

    // Select growth points on chains
    pair<Domain*, Domain*> growthpoint {select_new_growthpoint(selected_chain)};

    auto bound_domains {find_bound_domains(selected_chain)};
    unassign_domains(selected_chain);

    // Grow staple
    set_growthpoint_and_grow_staple(growthpoint, selected_chain);
    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    //HACK
    update_bias(+1);
    // Regrow staple in old conformation
    setup_for_regrow_old();

    // Select growth point from previously bound domains
    growthpoint = select_old_growthpoint(bound_domains);

    // Unassign and add to reversion list
    unassign_domains(selected_chain);

    // Grow staple
    set_growthpoint_and_grow_staple(growthpoint, selected_chain);

    // Revert modifier and test acceptance
    m_modifier = m_new_modifier;
    //HACK
    update_bias(+1);
    accepted = test_cb_acceptance();
    return accepted;
}

void CBStapleRegrowthMCMovetype::set_growthpoint_and_grow_staple(
        pair<Domain*, Domain*> growthpoint, vector<Domain*> selected_chain) {

    if (m_regrow_old) {
        set_old_growth_point(*growthpoint.first, *growthpoint.second);
    }
    else {
        double delta_e {set_growth_point(*growthpoint.first, *growthpoint.second)};
        m_bias *= exp(-delta_e);
    }
    if (not m_rejected) {
        grow_staple(growthpoint.first->m_d, selected_chain);
    }
}

void CBStapleRegrowthMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        select_and_set_config(domains, i);
        if (m_rejected) {
            break;
        }
    }
}

vector<double> CBStapleRegrowthMCMovetype::calc_bias(vector<double> bfactors,
        Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree,
        vector<Domain*>) {
    // Calculate rosenbluth weight
    double rosenbluth_i {0};
    for (auto bfactor: bfactors) {
        rosenbluth_i += bfactor;
    }
    vector<double> weights {};
    if (rosenbluth_i == 0) {
        
        // Deadend
        m_rejected = true;
    }
    else {
        m_bias *= rosenbluth_i;
        for (auto bfactor: bfactors) {
            weights.push_back(bfactor / rosenbluth_i);
        }
    }

    return weights;
}

bool CTCBScaffoldRegrowthMCMovetype::attempt_move() {
    bool accepted;

    m_regrow_old = false;
    vector<Domain*> scaffold_domains {select_scaffold_indices()};

    m_constraintpoints.calculate_constraintpoints(scaffold_domains);
    set<int> staples {m_constraintpoints.staples_to_be_regrown()};
    
    // Unassign staples except those directly and indirectly bound to external scaffold domains
    vector<Domain*> scaffold_domains_to_unassign(scaffold_domains.begin() + 1,
            scaffold_domains.end());
    unassign_domains(scaffold_domains_to_unassign);
    for (auto c_i: staples) {
        unassign_domains(m_origami_system.get_chain(c_i));
    }

    // Grow scaffold and staples
    if (m_constraintpoints.is_growthpoint(scaffold_domains[0])) {
        grow_staple_and_update_endpoints(scaffold_domains[0]);
        if (m_rejected) {
            accepted = false;
            return accepted;
        }
    }
    grow_chain(scaffold_domains);
    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    //HACK
    update_bias(+1);
    // Regrow in old conformation
    setup_for_regrow_old();
    m_constraintpoints.reset_active_endpoints();

    // Unassign staples except those directly and indirectly bound to external scaffold domains
    unassign_domains(scaffold_domains_to_unassign);
    for (auto c_i: staples) {
        unassign_domains(m_origami_system.get_chain(c_i));
    }

    // Grow scaffold and staples
    if (m_constraintpoints.is_growthpoint(scaffold_domains[0])) {
        grow_staple_and_update_endpoints(scaffold_domains[0]);
    }

    grow_chain(scaffold_domains);

    // Reset modifier and test acceptance
    m_modifier = 1;
    //HACK
    update_bias(+1);
    accepted = test_cb_acceptance();
    return accepted;
}

void CTCBScaffoldRegrowthMCMovetype::reset_internal() {
    CBMCMovetype::reset_internal();
    m_constraintpoints.reset_internal();
}

vector<Domain*> CTCBScaffoldRegrowthMCMovetype::select_scaffold_indices() {

    // Randomly select endpoints
    vector<Domain*> scaffold {m_origami_system.get_chain(m_origami_system.c_scaffold)};
    Domain* start_domain {scaffold[m_random_gens.uniform_int(0, scaffold.size() - 1)]};
    Domain* end_domain {scaffold[m_random_gens.uniform_int(0, scaffold.size() - 1)]};
    while (start_domain->m_d == end_domain->m_d) {
        end_domain = scaffold[m_random_gens.uniform_int(0, scaffold.size() - 1)];
    }
    
    vector<Domain*> domains {};

    // Cyclic domains
    if (m_origami_system.m_cyclic) {

        // Select direction of regrowth MESSY
        m_dir = m_random_gens.uniform_int(0, 1);
        if (m_dir == 0) {
            m_dir = -1;
        }
        Domain* endpoint_domain {(*end_domain) + m_dir};
        int d_i {start_domain->m_d};
        Domain* cur_domain {start_domain};
        while (d_i != end_domain->m_d) {
            domains.push_back(cur_domain);
            cur_domain = (*cur_domain) + m_dir;
            d_i = cur_domain->m_d;
        }
        domains.push_back(end_domain);
        m_constraintpoints.add_active_endpoint(endpoint_domain, endpoint_domain->m_pos);
    }

    // Linear domains
    else {

        // Find direction of regrowth
        if (end_domain->m_d > start_domain->m_d) {
            m_dir = 1;
        }
        else {
            m_dir = -1;
        }
        Domain* endpoint_domain {(*end_domain) + m_dir};
        for (int d_i {start_domain->m_d}; d_i != end_domain->m_d + m_dir; d_i += m_dir) {
            Domain* cur_domain {scaffold[d_i]};
            domains.push_back(cur_domain);
        }

        // If end domain is end of chain, no endpoint
        if (endpoint_domain != nullptr) {
            m_constraintpoints.add_active_endpoint(endpoint_domain, endpoint_domain->m_pos);
        }
    }

    return domains;
}

void CTCBScaffoldRegrowthMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        select_and_set_config(domains, i);
        if (m_rejected) {
            break;
        }
        Domain* domain {domains[i]};
        m_constraintpoints.update_endpoints(domain);
    
        // Grow staple if growth point
        if (m_constraintpoints.is_growthpoint(domain)) {
            grow_staple_and_update_endpoints(domain);
            if (m_rejected) {
                break;
            }
        }
    }
}

void CTCBScaffoldRegrowthMCMovetype::grow_staple_and_update_endpoints(
        Domain* growth_d_old) {

    Domain* growth_d_new {m_constraintpoints.get_domain_to_grow(growth_d_old)};
    int c_i {growth_d_new->m_c};
    if (m_regrow_old) {
        set_old_growth_point(*growth_d_new, *growth_d_old);
    }
    else {
        double delta_e {set_growth_point(*growth_d_new, *growth_d_old)};
        m_bias *= exp(-delta_e);
    }
    if (not m_rejected) {
        m_constraintpoints.update_endpoints(growth_d_new);
        vector<Domain*> staple {m_origami_system.get_chain(c_i)};
        grow_staple(growth_d_new->m_d, staple);
    }
}

vector<double> CTCBScaffoldRegrowthMCMovetype::calc_bias(
        vector<double> bfactors,
        Domain* domain,
        vector<pair<VectorThree, VectorThree>>& configs,
        VectorThree p_prev,
        vector<Domain*> domains) {


    vector<long double> weights(bfactors.begin(), bfactors.end());
    for (size_t i {0}; i != configs.size(); i++) {
        VectorThree cur_pos {configs[i].first};

        // Set weights of positions involving non-self binding to 0 unless endpoint reached
        Occupancy pos_occ {m_origami_system.position_occupancy(cur_pos)};
        if (pos_occ == Occupancy::unbound) {
            Domain* occ_domain {m_origami_system.unbound_domain_at(cur_pos)};
            bool binding_same_chain {occ_domain->m_c == domain->m_c};
            bool endpoint {m_constraintpoints.endpoint_reached(domain, cur_pos)};
            if (not (binding_same_chain or endpoint)) {
                weights[i] = 0;
            }
        }

        // Bias weights with number of walks
        weights[i] *= m_constraintpoints.calc_num_walks_prod(domain, cur_pos,
                domains, m_dir);
    }

    // Calculate number of walks for previous position
    long double prod_num_walks {m_constraintpoints.calc_num_walks_prod(domain,
            p_prev, domains, m_dir, 1)};

    // Modified Rosenbluth
    long double weights_sum {0};
    for (auto weight: weights) {
        weights_sum += weight;
    }

    // Check for deadend
    vector<double> norm_weights {};
    if (weights_sum == 0) {
        m_rejected = true;
    }
    else {
        long double bias {weights_sum / prod_num_walks};
        m_bias *= bias;

        // Normalize
        for (size_t i {0}; i != weights.size(); i++) {
            double norm_weight = static_cast<double>(weights[i] / weights_sum);
            norm_weights.push_back(norm_weight);
        }
    }

    return norm_weights;
}
