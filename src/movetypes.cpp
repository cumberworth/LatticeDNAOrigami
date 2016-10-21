// movetypes.cpp

#include <algorithm>
#include <random>
#include <set>

#include "utility.h"
#include "random_gens.h"
#include "movetypes.h"
#include "ideal_random_walk.h"

using std::fmin;
using std::set;
using std::find;

using namespace Movetypes;
using namespace Utility;
using namespace RandomGen;
using namespace IdealRandomWalk;

void MCMovetype::reset_origami() {
    
    // Have to deal with added domains, modified domains, and deleted domains

    // Unassign assigned domains
    for (auto c_i_d_i: m_assigned_domains) {
        int c_i {c_i_d_i.first};
        int d_i {c_i_d_i.second};
        pair<int, int> key {c_i, d_i};
        Domain* domain {m_origami_system.get_domain(c_i, d_i)};
        m_origami_system.unassign_domain(*domain);
    }

    // Delete chains that were added
    for (auto c_i: m_added_chains) {
        m_origami_system.delete_chain(c_i);

        // Roll back index numbering to prevent excessivly large indices
        m_origami_system.m_current_c_i -= 1;
    }

    // Revert modified domains to previous positions
    for (auto c_i_d_i: m_modified_domains) {
        int c_i {c_i_d_i.first};
        int d_i {c_i_d_i.second};
        pair<int, int> key {c_i, d_i};
        VectorThree pos {m_prev_pos[key]};
        VectorThree ore {m_prev_ore[key]};
        Domain* domain {m_origami_system.get_domain(c_i, d_i)};
        m_origami_system.set_checked_domain_config(*domain, pos, ore);
    }
    m_origami_system.m_constraints_violated = false;
}

Domain* MCMovetype::select_random_domain() {
    int d_i_index {m_random_gens.uniform_int(0, m_origami_system.num_domains() - 1)};
    int counter {0};
    for (auto chain: m_origami_system.get_chains()) {
        for (auto domain: chain) {
            if (counter == d_i_index) {
                return domain;
            }
            else {
                counter++;
            }
        }
    }
    return nullptr;
}

int MCMovetype::select_random_staple_identity() {
    return m_random_gens.uniform_int(1, m_origami_system.m_identities.size() - 1);
}

int MCMovetype::select_random_staple_of_identity(int c_i_ident) {
    int staples {m_origami_system.num_staples_of_ident(c_i_ident)};
    if (staples == 0) {
        m_rejected = true;
        return -1;
    }
    int staple_i {m_random_gens.uniform_int(0, staples - 1)};
    int c_i {m_origami_system.staples_of_ident(c_i_ident)[staple_i]};

    return c_i;
}

VectorThree MCMovetype::select_random_position(VectorThree p_prev) {
    VectorThree vec {vectors[m_random_gens.uniform_int(0, 5)]};
    return p_prev + vec;
}

VectorThree MCMovetype::select_random_orientation() {
    VectorThree vec {vectors[m_random_gens.uniform_int(0, 5)]};
    return vec;
}

bool MCMovetype::test_acceptance(double p_ratio) {
    double p_accept = fmin(1, p_ratio) * m_modifier;
    bool accept;
    if (p_accept == 1) {
        accept = true;
    }
    else {
        if (p_accept > m_random_gens.uniform_real()) {
            accept = true;
        }
        else {
            accept = false;
        }
    }

    return accept;
}

bool MCMovetype::staple_is_connector(vector<Domain*> staple) {
    for (auto domain: staple) {
        if (domain->m_state != Occupancy::unbound) {
            Domain* bound_domain {domain->m_bound_domain};
            if (bound_domain->m_c == m_origami_system.c_scaffold) {
                continue;
            }

            // Do I really need to clear this?
            set<int> participating_chains {domain->m_c};
            if (not scan_for_scaffold_domain(bound_domain, participating_chains)) {
                return true;
            }
        }
    }
    return false;
}

bool MCMovetype::scan_for_scaffold_domain(
        Domain* domain,
        set<int>& participating_chains) {

    bool bound_to_scaffold {false};
    int c_i {domain->m_c};
    participating_chains.insert(c_i);
    vector<Domain*> staple {m_origami_system.get_chain(c_i)};
    for (auto cur_domain: staple) {
        if (cur_domain == domain) {
            continue;
        }
        Domain* bound_domain {cur_domain->m_bound_domain};
        if (bound_domain != nullptr) {

            // Skip if bound to self
            if (bound_domain->m_c == domain->m_c) {
                continue;
            }

            // Check if bound to scaffold
            bool domain_on_scaffold {bound_domain->m_c == m_origami_system.c_scaffold};
            if (domain_on_scaffold) {
                return true;
            }

            // Check if bound domain already in progress
            if (participating_chains.count(bound_domain->m_c) > 0) {
                continue;
            }
            else {
                bound_to_scaffold = scan_for_scaffold_domain(bound_domain,
                        participating_chains);
            }
            if (bound_to_scaffold) {
                return true;
            }
        }
    }
    return false;
}

double RegrowthMCMovetype::set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old) {
    // Set growth point with new orientation
    VectorThree o_new = select_random_orientation();
    double delta_e {0};
    delta_e += m_origami_system.set_domain_config(growth_domain_new,
            growth_domain_old.m_pos, o_new);
    if (m_origami_system.m_constraints_violated) {
        m_rejected = true;
    }
    else {
        pair<int, int> key {growth_domain_new.m_c, growth_domain_new.m_d};
        m_assigned_domains.push_back(key);
    }

    return delta_e;
}

void RegrowthMCMovetype::grow_staple(int d_i_index, vector<Domain*> selected_chain) {
    // Grow staple in both directions out from growth point
    // The indices passed to grow chain should include the growth point
    int old_num_bound_domains {m_origami_system.num_fully_bound_domain_pairs() -
            m_origami_system.num_self_bound_domain_pairs()};

    // Grow in three prime direction (staple domains increase in 3' direction)
    auto first_iter3 {selected_chain.begin() + d_i_index};
    auto last_iter3 {selected_chain.end()};
    vector<Domain*> domains_three_prime {first_iter3, last_iter3};
    grow_chain(domains_three_prime);
    if (m_rejected) {
        return;
    }

    // Grow in five prime direction
    auto first_iter5 {selected_chain.begin()};
    auto last_iter5 {selected_chain.begin() + d_i_index + 1};
    vector<Domain*> domains_five_prime {first_iter5, last_iter5};
    std::reverse(domains_five_prime.begin(), domains_five_prime.end());
    grow_chain(domains_five_prime);
    if (m_rejected) {
        return;
    }

    // Overcount correction
    int overcount {m_origami_system.num_fully_bound_domain_pairs() -
            m_origami_system.num_self_bound_domain_pairs() -
            old_num_bound_domains};
    m_modifier /= (overcount + 1);
}

pair<Domain*, Domain*> RegrowthMCMovetype::select_new_growthpoint(
        vector<Domain*> selected_chain) {

    int growth_d_new {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_new {selected_chain[growth_d_new]};
    Domain* growth_domain_old {select_random_domain()};
    while (growth_domain_old->m_c == growth_domain_new->m_c) {
        growth_domain_old = select_random_domain();
    }
    return {growth_domain_new, growth_domain_old};
}

bool OrientationRotationMCMovetype::attempt_move() {
    bool accept;
    
    // Select random chain, domain, and orientation
    Domain* domain {select_random_domain()};
    VectorThree o_new {select_random_orientation()};

    if (domain->m_state == Occupancy::bound) {
        VectorThree o_old {domain->m_ore};
        Domain* bound_domain {domain->m_bound_domain};
        m_origami_system.unassign_domain(*bound_domain);
        m_origami_system.set_domain_orientation(*domain, o_new);
        VectorThree pos {domain->m_pos};
        m_origami_system.set_domain_config(*bound_domain, pos, -o_new);
        if (m_origami_system.m_constraints_violated) {
            accept = false;
            m_origami_system.unassign_domain(*bound_domain);
            m_origami_system.set_domain_orientation(*domain, o_old);
            m_origami_system.set_checked_domain_config(*bound_domain, pos, -o_old);
        }
        else {
            accept = true;
        }
    }
    else {
        m_origami_system.set_domain_orientation(*domain, o_new);
        accept = true;
    }


    return accept;
}

void MetMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        Domain* domain {domains[i]};
        Domain* prev_domain {domains[i - 1]};
        VectorThree new_p {select_random_position(prev_domain->m_pos)};
        VectorThree new_o {select_random_orientation()};
        m_delta_e += m_origami_system.set_domain_config(*domain, new_p, new_o);
        if (m_origami_system.m_constraints_violated) {
            m_rejected = true;
            break;
        }
        else {
            pair<int, int> key {domain->m_c, domain->m_d};
            m_assigned_domains.push_back(key);
        }
    }
}

void MetMCMovetype::unassign_domains(vector<Domain*> domains) {
    for (auto domain: domains) {
        pair<int, int> key {domain->m_c, domain->m_d};
        m_prev_pos[key] = domain->m_pos;
        m_prev_ore[key] = domain->m_ore;
        m_modified_domains.push_back(key);
        m_delta_e += m_origami_system.unassign_domain(*domain);
    }
}

bool MetStapleExchangeMCMovetype::attempt_move() {

    // Select inertion or deletion with equal frequency
    bool accept;
    if (m_random_gens.uniform_real() < 0.5) {
        accept = insert_staple();
    }
    else {
        accept = delete_staple();
    }
    
    return accept;
}

bool MetStapleExchangeMCMovetype::staple_insertion_accepted(int c_i_ident) {
    double boltz_factor {exp(-m_delta_e)};
    int Ni_new {m_origami_system.num_staples_of_ident(c_i_ident)};

    // Correct for extra states from additional staple domains
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
    int extra_states {static_cast<int>(pow(6, extra_df))};
    double ratio {extra_states / Ni_new * boltz_factor};

    // Correct for insertion into subset of volume
    m_modifier *= m_insertion_sites / m_origami_system.m_volume;

    // Correct for considering only 1 of staple length ways insertion could occur
    m_modifier *= staple_length;

    return test_acceptance(ratio);
}

bool MetStapleExchangeMCMovetype::staple_deletion_accepted(int c_i_ident) {
    double boltz_factor {exp(-m_delta_e)};
    int Ni_new {m_origami_system.num_staples()};

    // Correct for extra states from additional staple domains
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
    int extra_states {static_cast<int>(pow(6, extra_df))};

    double ratio {(Ni_new + 1) / extra_states * boltz_factor};

    return test_acceptance(ratio);
}

bool MetStapleExchangeMCMovetype::insert_staple() {
    bool accepted;

    // Select and add chain of random identity
    int c_i_ident {select_random_staple_identity()};
    int c_i {m_origami_system.add_chain(c_i_ident)};
    m_added_chains.push_back(c_i);

    // Assume that add_chain always adds to end of m_domains
    vector<Domain*> selected_chain {m_origami_system.get_last_chain()};

    // Select growth point on added chain
    int d_i_index {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_new {selected_chain[d_i_index]};

    // Select and set growth point on current system
    select_and_set_growth_point(growth_domain_new);
    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    grow_staple(d_i_index, selected_chain);
    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    accepted = staple_insertion_accepted(c_i_ident);
    return accepted;
}

bool MetStapleExchangeMCMovetype::delete_staple() {
    bool accepted;

    cout << "BROKEN\n";
    // Select random identity and staple of that type
    int c_i_ident {select_random_staple_identity()};
    int c_i {select_random_staple_of_identity(c_i_ident)};
    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    // Add stuff to reversions lists and unassign domsins
    vector<Domain*> staple {m_origami_system.get_chain(c_i)};
    int chain_length {static_cast<int>(staple.size())};
    for (int d_i {0}; d_i != chain_length; d_i++) {
        Domain* domain {staple[d_i]};
        pair<int, int> key {c_i, d_i};
        m_prev_pos[key] = domain->m_pos;
        m_prev_ore[key] = domain->m_ore;
        m_delta_e += m_origami_system.unassign_domain(*domain);
    }

    // Delete chain and test acceptance
    m_origami_system.delete_chain(c_i);

    return staple_deletion_accepted(c_i_ident);
}

void MetStapleExchangeMCMovetype::select_and_set_growth_point(Domain* growth_domain_new) {

    // Select a different chain, excluding newest addition which is last index
    int c_i_index {m_random_gens.uniform_int(0, m_origami_system.num_staples() - 1)};
    vector<Domain*> selected_chain {m_origami_system.get_chains()[c_i_index]};

    // Select a domain of the chain
    int d_i  {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_old {selected_chain[d_i]};

    m_delta_e += set_growth_point(*growth_domain_new, *growth_domain_old);
}

bool MetStapleRegrowthMCMovetype::attempt_move() {
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

    unassign_domains(selected_chain);

    // Grow staple
    m_delta_e += set_growth_point(*growthpoint.first, *growthpoint.second);
    if (m_rejected) {
        accepted = false;
        return accepted;
    }
    grow_staple(growthpoint.first->m_d, selected_chain);
    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    double boltz_factor {exp(-m_delta_e)};
    accepted = test_acceptance(boltz_factor);
    return accepted;
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
                bool domains_complementary {m_origami_system.check_domains_complementary(
                        domain, *unbound_domain)};
                VectorThree o_new;
                double bfactor;
                if (domains_complementary) {
                    o_new = -unbound_domain->m_ore;
                    bfactor = 1;
                }
                else {
                    o_new = {0, 0, 0};
                    bfactor = 6;
                }
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
                bfactors.push_back(6);
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

double CBMCMovetype::set_old_growth_point(Domain& growth_domain_new, Domain& growth_domain_old) {
    pair<int, int> key {growth_domain_new.m_c, growth_domain_new.m_d};
    VectorThree o_old {m_old_ore[key]};
    double delta_e {0};
    delta_e += m_origami_system.set_checked_domain_config(growth_domain_new,
            growth_domain_old.m_pos, o_old);
    m_assigned_domains.push_back(key);

    return delta_e;
}

bool CBMCMovetype::test_cb_acceptance() {
    double ratio {m_new_bias / m_bias};
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

bool CBStapleExchangeMCMovetype::attempt_move() {
    // Select inertion or deletion with equal frequency
    bool accept;
    if (m_random_gens.uniform_real() < 0.5) {
        accept = insert_staple();
    }
    else {
        accept = delete_staple();
    }
    
    return accept;
}

double CBStapleExchangeMCMovetype::calc_staple_insertion_acc_ratio(int c_i_ident) {
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    m_bias /= pow(6, staple_length);
    int Ni_new {m_origami_system.num_staples_of_ident(c_i_ident)};

    // Correct for extra states from additional staple domains
    int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
    double extra_states {pow(6, extra_df)};
    double ratio {extra_states / static_cast<double>(Ni_new) * m_bias};

    // Correct for insertion into subset of volume
    m_modifier *= m_insertion_sites / m_origami_system.m_volume;

    // Correct for considering only 1 of staple length ways insertion could occur
    m_modifier *= staple_length;

    return ratio;
}

double CBStapleExchangeMCMovetype::calc_staple_deletion_acc_ratio(int c_i_ident) {
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    m_bias /= pow(6, staple_length);
    int Ni {m_origami_system.num_staples_of_ident(c_i_ident)};

    // Correct for extra states from additional staple domains
    int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
    double extra_states {pow(6, extra_df)};

    double ratio {static_cast<double>(Ni) / extra_states / m_bias};

    return ratio;
}

vector<double> CBStapleExchangeMCMovetype::calc_bias(vector<double> weights,
        Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree,
        vector<Domain*>) {
    // Calculate rosenbluth weight
    double rosenbluth_i {0};
    for (auto weight: weights) {
        rosenbluth_i += weight;
    }
    if (rosenbluth_i == 0) {
        
        // Deadend
        m_rejected = true;
    }
    else {
        m_bias *= rosenbluth_i;
        for (size_t i {0}; i != weights.size(); i++) {
            weights[i] /= rosenbluth_i;
        }
    }

    return weights;
}

bool CBStapleExchangeMCMovetype::insert_staple() {
    bool accepted;

    // Select and add chain of random identity
    int c_i_ident {select_random_staple_identity()};
    int c_i {m_origami_system.add_chain(c_i_ident)};
    m_added_chains.push_back(c_i);

    // Assume that add_chain always adds to end of m_domains
    vector<Domain*> selected_chain {m_origami_system.get_last_chain()};

    // Select growth point on added chain
    int d_i {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_new {selected_chain[d_i]};

    // Select and set growth point on current system
    select_and_set_growth_point(growth_domain_new);
    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    grow_staple(d_i, selected_chain);
    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    double ratio {calc_staple_insertion_acc_ratio(c_i_ident)};
    accepted = test_acceptance(ratio);
    return accepted;
}

bool CBStapleExchangeMCMovetype::delete_staple() {
    bool accepted;

    // Select random identity and staple of that type
    int c_i_ident {select_random_staple_identity()};
    int c_i {select_random_staple_of_identity(c_i_ident)};
    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    // Reject if staple is connector
    vector<Domain*> staple {m_origami_system.get_chain(c_i)};
    if (staple_is_connector(staple)) {
        accepted = false;
        return accepted;
    }

    // Select growth point
    auto bound_domains {find_bound_domains(staple)};
    auto growthpoint {select_old_growthpoint(bound_domains)};
    
    // Add stuff to reversions lists and unassign domains
    unassign_for_regrowth(staple);

    // Regrow chain to calculate Rosenbluth weight
    m_regrow_old = true;
    double delta_e {set_old_growth_point(*growthpoint.first, *growthpoint.second)};
    m_bias *= 6 * exp(-delta_e);
    grow_staple(growthpoint.first->m_d, staple);

    // Remove overcounting correction
    m_modifier = 1;
    double ratio {calc_staple_deletion_acc_ratio(c_i_ident)};
    accepted = test_acceptance(ratio);

    if (accepted) {
        unassign_and_delete_staple(c_i, staple);
    }

    // No need for reseting with this movetype
    m_assigned_domains.clear();

    return accepted;
}

void CBStapleExchangeMCMovetype::unassign_and_delete_staple(
        int c_i,
        vector<Domain*> staple) {
    for (auto domain: staple) {
        m_origami_system.unassign_domain(*domain);
    }
    m_origami_system.delete_chain(c_i);
}

void CBStapleExchangeMCMovetype::unassign_for_regrowth(vector<Domain*> domains) {
    // No need for reversion list
    for (auto domain: domains) {
        pair<int, int> key {domain->m_c, domain->m_d};
        m_old_pos[key] = domain->m_pos;
        m_old_ore[key] = domain->m_ore;

        // Would be double counting to include delta_e in bias
        m_origami_system.unassign_domain(*domain);
    }
}

void CBStapleExchangeMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        select_and_set_config(domains, i);
        if (m_rejected) {
            break;
        }
    }
}

void CBStapleExchangeMCMovetype::select_and_set_growth_point(Domain* growth_domain_new) {

    // Select a different chain, excluding newest addition which is last index
    int c_i_index {m_random_gens.uniform_int(0, m_origami_system.num_staples() - 1)};
    vector<Domain*> selected_chain {m_origami_system.get_chains()[c_i_index]};

    // Select a domain of the chain
    int d_i  {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_old {selected_chain[d_i]};

    double delta_e {set_growth_point(*growth_domain_new, *growth_domain_old)};
    if (m_origami_system.m_constraints_violated) {
        m_rejected = true;
    }
    m_bias *= 6 * exp(-delta_e);
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
    accepted = test_cb_acceptance();
    return accepted;
}

void CBStapleRegrowthMCMovetype::set_growthpoint_and_grow_staple(
        pair<Domain*, Domain*> growthpoint, vector<Domain*> selected_chain) {

    double delta_e {0};
    if (m_regrow_old) {
        delta_e = set_old_growth_point(*growthpoint.first, *growthpoint.second);
    }
    else {
        delta_e = set_growth_point(*growthpoint.first, *growthpoint.second);
    }
    if (not m_rejected) {
        m_bias *= exp(-delta_e);
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

Constraintpoints::Constraintpoints(
        OrigamiSystem& origami_system,
        IdealRandomWalks& ideal_random_walks) :
        m_origami_system {origami_system},
        m_ideal_random_walks {ideal_random_walks} {
}

void Constraintpoints::calculate_constraintpoints(
        vector<Domain*> scaffold_domains) {
    find_staples_growthpoints_endpoints(scaffold_domains);

    // Save initial active endpionts and remaining steps for regrowing old
    m_initial_active_endpoints = m_active_endpoints;
}

bool Constraintpoints::is_growthpoint(Domain* domain) {
    return (m_growthpoints.count(domain) > 0);
}

void Constraintpoints::add_active_endpoint(Domain* domain, VectorThree pos) {
    m_active_endpoints[domain->m_c].push_back({domain->m_d, pos});
}

void Constraintpoints::reset_active_endpoints() {
    m_active_endpoints = m_initial_active_endpoints;
}

void Constraintpoints::remove_active_endpoint(Domain* domain) {

    // Remove endpoint if reached
    int c_i {domain->m_c};
    vector<int> endpoints_to_erase {};
    for (size_t j {0}; j != m_active_endpoints[c_i].size(); j++) {
        pair<int, VectorThree> endpoint {m_active_endpoints[c_i][j]};
        if (endpoint.first == domain->m_d) {
            endpoints_to_erase.push_back(j);
        }
    }

    // Removing in reverse means the indices are not changing as I do it
    std::reverse(endpoints_to_erase.begin(), endpoints_to_erase.end());
    for (auto j: endpoints_to_erase) {
        m_active_endpoints[c_i].erase(m_active_endpoints[c_i].begin() + j);
    }
}

void Constraintpoints::update_endpoints(Domain* domain) {
    remove_active_endpoint(domain);

    // Copy inactive endpoints to active endpoints
    bool inactive_endpoint_present {m_inactive_endpoints.count(domain) > 0};
    if (inactive_endpoint_present) {
        Domain* endpoint_domain {m_inactive_endpoints[domain]};
        add_active_endpoint(endpoint_domain, domain->m_pos);
        // pretty sure it's fine not to erase the inverse inactive endpoint
    }
}


Domain* Constraintpoints::get_domain_to_grow(Domain* domain) {
    return m_growthpoints[domain];
}

bool Constraintpoints::endpoint_reached(Domain* domain, VectorThree pos) {
    bool reached {false};
    for (auto endpoint: m_active_endpoints[domain->m_c]) {
        bool endpoint_domain_reached {domain->m_d == endpoint.first};
        bool endpoint_pos_reached {pos == endpoint.second};
        if (endpoint_domain_reached and endpoint_pos_reached) {
            reached = true;
            break;
        }
    }
    return reached;
}

double Constraintpoints::calc_num_walks_prod(
        Domain* domain,
        VectorThree pos,
        vector<Domain*> domains,
        // For calculating previous position endpoints
        int step_offset) {

    double prod_num_walks {1};
    for (auto endpoint: m_active_endpoints[domain->m_c]) {

        // Check if endpoint in current stretch of domains being regrown (
        // staples are grown out from there growth pointion two seperate
        // directions
        int endpoint_d_i {endpoint.first};
        int first_d_i {domains.front()->m_d};
        int last_d_i {domains.back()->m_d};
        if (domain->m_c == m_origami_system.c_scaffold or
                (endpoint_d_i >= first_d_i and endpoint_d_i <= last_d_i) or
                (endpoint_d_i >= last_d_i and endpoint_d_i <= first_d_i)) {
            int steps {abs(endpoint_d_i - domain->m_d) + step_offset};
            VectorThree endpoint_p {endpoint.second};
            prod_num_walks *= m_ideal_random_walks.num_walks(pos, endpoint_p,
                    steps);
        }
    }
    return prod_num_walks;
}

void Constraintpoints::find_staples_growthpoints_endpoints(
        vector<Domain*> scaffold_domains) {
    
    // Chain ids of staples that have been checked
    set<int> checked_staples {};
    for (auto domain: scaffold_domains) {

        // Domains bound in network extending from current scaffold domain
        set<int> participating_chains {m_origami_system.c_scaffold};
        vector<pair<Domain*, Domain*>> potential_growthpoints {};
        bool externally_bound {false};
        bool bound {domain->m_bound_domain != nullptr};

        // Skip if not bound, bound to self or already checked
        if (not bound or (domain->m_bound_domain->m_c == domain->m_c) or
                staple_already_checked(domain, checked_staples)) {
            continue;
        }
        else {
            Domain* bound_domain {domain->m_bound_domain};
            potential_growthpoints.push_back({domain, bound_domain});

            // Finds network and updates given vectors and sets
            scan_staple_topology(bound_domain, participating_chains,
                    potential_growthpoints, scaffold_domains, externally_bound);

            if (externally_bound) {

                // Network will not be regrown
                // Convert inactive endpoints on scaffold to active
                add_active_endpoints_on_scaffold(participating_chains,
                        potential_growthpoints);
            }
            else {

                // Will need to grow network
                add_growthpoints(potential_growthpoints);
                add_regrowth_staples(participating_chains);
            }

            // Add to checked staples
            for (auto c_i: participating_chains) {
                checked_staples.insert(c_i);
            }
        }
    }

    // Remove endpoints on first scaffold domain
    remove_active_endpoint(scaffold_domains[0]);
}

bool Constraintpoints::staple_already_checked(
        Domain* domain,
        set<int> checked_staples) {
    Domain* bound_domain {domain->m_bound_domain};
    bool staple_checked {find(checked_staples.begin(), checked_staples.end(),
            bound_domain->m_c) != checked_staples.end()};
    return staple_checked;
}

void Constraintpoints::add_growthpoints(
        vector<pair<Domain*, Domain*>> potential_growthpoints) {

    for (auto growthpoint: potential_growthpoints) {
        m_growthpoints[growthpoint.first] = growthpoint.second;
    }
}

void Constraintpoints::add_regrowth_staples(set<int> participating_chains) {
    for (auto c_i: participating_chains) {
        if (c_i == 0) {
            continue;
        }
        m_regrowth_staples.insert(c_i);
    }
}

void Constraintpoints::add_active_endpoints_on_scaffold(
        set<int> participating_chains,
        vector<pair<Domain*, Domain*>> potential_growthpoints) {

    // Extract those from inactive endpoints
    for (auto c_i: participating_chains) {
        for (auto domain: m_origami_system.get_chain(c_i)) {
            if (m_inactive_endpoints.count(domain) > 0) {
                Domain* endpoint_domain {m_inactive_endpoints[domain]};
                if (endpoint_domain->m_c == 0) {
                    add_active_endpoint(endpoint_domain, endpoint_domain->m_pos);
                }
            // else delete endpoint (but unnecessary so don't bother)
            }
        }
    }

    // Extract those from potential growthpoints
    for (auto growth_pair: potential_growthpoints) {
        Domain* growth_domain {growth_pair.first};
        if (growth_domain->m_c == 0) {
            add_active_endpoint(growth_domain, growth_domain->m_pos);
        }
    }
}

void Constraintpoints::scan_staple_topology(
        Domain* growth_domain,
        set<int>& participating_chains,
        vector<pair<Domain*, Domain*>>& potential_growthpoints,
        vector<Domain*>& scaffold_domains,
        bool& externally_bound) {

    int c_i {growth_domain->m_c};
    participating_chains.insert(c_i);
    vector<Domain*> staple {m_origami_system.get_chain(c_i)};
    for (auto domain: staple) {
        if (growth_domain == domain) {
            continue;
        }
        Domain* bound_domain {domain->m_bound_domain};
        bool domain_bound {bound_domain != nullptr};

        // Skip if not bound or bound to self
        if (not domain_bound or (bound_domain->m_c == c_i)) {
            continue;
        }
        else {
            int bound_c_i {bound_domain->m_c};
            bool in_network {participating_chains.count(bound_c_i) > 0};
            if (in_network) {

                // Add domain to endpoints
                bool domain_in_range {find(scaffold_domains.begin(), scaffold_domains.end(),
                        bound_domain) != scaffold_domains.end()};
                if (not domain_in_range) {
                    externally_bound = true;
                }
                else {
                    m_inactive_endpoints[domain] = bound_domain;
                    m_inactive_endpoints[bound_domain] = domain;
                }
            }
            else {

                // Add domain to potential growthpoints and scan the chain
                potential_growthpoints.push_back({domain, bound_domain});
                scan_staple_topology(bound_domain, participating_chains,
                        potential_growthpoints, scaffold_domains, externally_bound);
            }
        }
    }
}

bool CTCBScaffoldRegrowthMCMovetype::attempt_move() {
    bool accepted;

    m_regrow_old = false;
    //DEBUG
    //vector<Domain*> scaffold_domains {select_scaffold_indices()};
    vector<Domain*> scaffold_domains {m_origami_system.get_chain(0)};

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
    accepted = test_cb_acceptance();
    return accepted;
}

vector<Domain*> CTCBScaffoldRegrowthMCMovetype::select_scaffold_indices() {

    // Randomly select endpoints
    vector<Domain*> scaffold {m_origami_system.get_chain(m_origami_system.c_scaffold)};
    Domain* start_domain {scaffold[m_random_gens.uniform_int(0, scaffold.size() - 1)]};
    Domain* end_domain {scaffold[m_random_gens.uniform_int(0, scaffold.size() - 1)]};
    while (start_domain->m_d == end_domain->m_d) {
        end_domain = scaffold[m_random_gens.uniform_int(0, scaffold.size() - 1)];
    }
    
    // Find direction of regrowth
    int dir;
    if (end_domain->m_d > start_domain->m_d) {
        dir = 1;
    }
    else {
        dir = -1;
    }
    Domain* endpoint_domain {(*end_domain) + dir};

    vector<Domain*> domains {};

    // Cyclic domains
    if (m_origami_system.m_cyclic) {
        int d_i {start_domain->m_d};
        Domain* cur_domain {start_domain};
        while (d_i != end_domain->m_d) {
            domains.push_back(cur_domain);
            cur_domain = cur_domain + 1;
            d_i = cur_domain->m_d;
        }
        domains.push_back(end_domain);
        m_constraintpoints.add_active_endpoint(endpoint_domain, endpoint_domain->m_pos);
    }

    // Linear domains
    else {
        for (int d_i {start_domain->m_d}; d_i != end_domain->m_d + dir; d_i += dir) {
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
    double delta_e;
    if (m_regrow_old) {
        delta_e = set_old_growth_point(*growth_d_new, *growth_d_old);
    }
    else {
        delta_e = set_growth_point(*growth_d_new, *growth_d_old);
    }
    if (not m_rejected) {
        m_constraintpoints.update_endpoints(growth_d_new);
        m_bias *= exp(-delta_e);
        vector<Domain*> staple {m_origami_system.get_chain(c_i)};
        grow_staple(growth_d_new->m_d, staple);
    }
}

vector<double> CTCBScaffoldRegrowthMCMovetype::calc_bias(
        vector<double> weights,
        Domain* domain,
        vector<pair<VectorThree, VectorThree>>& configs,
        VectorThree p_prev,
        vector<Domain*> domains) {

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
        //DEBUG
        //weights[i] *= m_constraintpoints.calc_num_walks_prod(domain, cur_pos,
        //        domains);
    }

    // Calculate number of walks for previous position
    //DEBUG
    //double prod_num_walks {m_constraintpoints.calc_num_walks_prod(domain,
     //       p_prev, domains, 1)};

    // Modified Rosenbluth
    double weights_sum {0};
    for (auto weight: weights) {
        weights_sum += weight;
    }

    // Check for deadend
    if (weights_sum == 0) {
        m_rejected = true;
    }
    else {
        //DEBUG
        //double bias {weights_sum / prod_num_walks};
        double bias {weights_sum};
        m_bias *= bias;

        // Normalize
        for (size_t i {0}; i != weights.size(); i++) {
            weights[i] /= weights_sum;
        }
    }

    return weights;
}

template<typename T>
unique_ptr<MCMovetype> Movetypes::movetype_constructor(OrigamiSystem& origami_system, RandomGens& random_gens, IdealRandomWalks& ideal_random_walks) {
    return unique_ptr<MCMovetype> {new T {origami_system, random_gens, ideal_random_walks}};
}
