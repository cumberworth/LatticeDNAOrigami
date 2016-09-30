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
        int c_i_index {index(m_origami_system.m_chain_indices, c_i)};
        pair<int, int> key {c_i, d_i};
        Domain* domain {m_origami_system.m_domains[c_i_index][d_i]};
        m_origami_system.unassign_domain(*domain);
    }

    // Delete chains that were added
    for (auto c_i: m_added_chains) {
        m_origami_system.delete_chain(c_i);

        // Roll back index numbering to prevent excessivly large indices
        m_origami_system.m_current_c_i -= 1;
    }

    // Re-add chains that were deleted and revert domains to previous positions
    for (auto c_index_identity: m_deleted_chains) {
        int c_i {c_index_identity.first};
        m_origami_system.add_chain(c_index_identity.second, c_i);
        for (auto d_i: m_deleted_domains[c_i]) {
            int c_i_index {index(m_origami_system.m_chain_indices, c_i)};
            pair<int, int> key {c_i, d_i};
            VectorThree pos {m_prev_pos[key]};
            VectorThree ore {m_prev_ore[key]};
            Domain* domain {m_origami_system.m_domains[c_i_index][d_i]};
            m_origami_system.set_checked_domain_config(*domain, pos, ore);
        }
    }

    // Revert modified domains to previous positions
    for (auto c_i_d_i: m_modified_domains) {
        int c_i {c_i_d_i.first};
        int d_i {c_i_d_i.second};
        int c_i_index {index(m_origami_system.m_chain_indices, c_i)};
        pair<int, int> key {c_i, d_i};
        VectorThree pos {m_prev_pos[key]};
        VectorThree ore {m_prev_ore[key]};
        Domain* domain {m_origami_system.m_domains[c_i_index][d_i]};
        m_origami_system.set_checked_domain_config(*domain, pos, ore);
    }
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

Domain* MCMovetype::select_random_domain() {
    int d_i_index {m_random_gens.uniform_int(0, m_origami_system.num_domains() - 1)};
    int counter {0};
    for (auto chain: m_origami_system.m_domains) {
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
    size_t staples {m_origami_system.m_identity_to_index[c_i_ident].size()};
    if (staples == 0) {
        throw MoveRejection {};
    }
    int staple_i {m_random_gens.uniform_int(0, staples - 1)};
    int c_i {m_origami_system.m_identity_to_index[c_i_ident][staple_i]};

    return c_i;
}

bool MCMovetype::staple_is_connector(vector<Domain*> staple) {
    for (auto cur_domain: staple) {
        if (cur_domain->m_state != Occupancy::unbound) {
            Domain* bound_domain {cur_domain->m_bound_domain};
            if (bound_domain->m_c == m_origami_system.c_scaffold) {
                continue;
            }

            // Do I really need to clear this?
            set<int> participating_chains {cur_domain->m_c};
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
    participating_chains.insert(domain->m_c);
    int c_i_index {index(m_origami_system.m_chain_indices, domain->m_c)};
    vector<Domain*> staple {m_origami_system.m_domains[c_i_index]};
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
    try {
        delta_e += m_origami_system.set_domain_config(growth_domain_new,
                growth_domain_old.m_pos, o_new);
        pair<int, int> key {growth_domain_new.m_c, growth_domain_new.m_d};
        m_assigned_domains.push_back(key);
    }
    catch (ConstraintViolation) {
        throw MoveRejection {};
    }

    return delta_e;
}

void RegrowthMCMovetype::grow_staple(int d_i_index, vector<Domain*> selected_chain) {
    // Grow staple in both directions out from growth point
    // The indices passed to grow chain should include the growth point
    int old_num_bound_domains {m_origami_system.num_domains()};

    // Grow in three prime direction (staple domains increase in 3' direction)
    auto first_iter3 {selected_chain.begin() + d_i_index};
    auto last_iter3 {selected_chain.end()};
    vector<Domain*> domains_three_prime {first_iter3, last_iter3};
    grow_chain(domains_three_prime);

    // Grow in five prime direction
    auto first_iter5 {selected_chain.begin()};
    auto last_iter5 {selected_chain.begin() + d_i_index + 1};
    vector<Domain*> domains_five_prime {first_iter5, last_iter5};
    std::reverse(domains_five_prime.begin(), domains_five_prime.end());
    grow_chain(domains_five_prime);

    // Overcount correction
    int overcount {old_num_bound_domains - m_origami_system.num_domains()};
    m_modifier /= (overcount + 1);
}

bool OrientationRotationMCMovetype::attempt_move() {
    bool accept;
    
    // Select random chain and domain
    Domain* domain {select_random_domain()};

    // Reject if in bound state
    if (domain->m_state == Occupancy::bound) {
        accept = false;
        return accept;
    }

    // All other cases accept
    else {
        accept = true;
    }

    // Select random orientation and update
    VectorThree o_new {select_random_orientation()};
    m_origami_system.set_domain_orientation(*domain, o_new);

    return accept;
}

void MetMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        Domain* domain {domains[i]};
        Domain* prev_domain {domains[i - 1]};
        VectorThree new_p {select_random_position(prev_domain->m_pos)};
        VectorThree new_o {select_random_orientation()};
        try {
            m_origami_system.set_domain_config(*domain, new_p, new_o);
            pair<int, int> key {domain->m_c, domain->m_d};
            m_assigned_domains.push_back(key);
        }
        catch (ConstraintViolation) {
            throw MoveRejection {};
        }
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
    int Ni_new {m_origami_system.num_staples()};

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

    // Select and add chain of random identity
    int c_i_ident {select_random_staple_identity()};
    int c_i {m_origami_system.add_chain(c_i_ident)};
    m_added_chains.push_back(c_i);

    // Assume that add_chain always adds to end of m_domains
    vector<Domain*> selected_chain {m_origami_system.m_domains.back()};

    // Select growth point on added chain
    int d_i_index {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_new {selected_chain[d_i_index]};

    // Select and set growth point on current system
    select_and_set_growth_point(growth_domain_new);

    grow_staple(d_i_index, selected_chain);

    return staple_insertion_accepted(c_i_ident);
}

bool MetStapleExchangeMCMovetype::delete_staple() {

    // Select random identity and staple of that type
    int c_i_ident {select_random_staple_identity()};
    int c_i {select_random_staple_of_identity(c_i_ident)};

    // The index into m_domains
    int c_i_index {index(m_origami_system.m_chain_indices, c_i)};
    
    // Add stuff to reversions lists and unassign domsins
    int chain_length {static_cast<int>(m_origami_system.m_domains[c_i_index].size())};
    for (int d_i {0}; d_i != chain_length; d_i++) {
        Domain* domain {m_origami_system.m_domains[c_i_index][d_i]};
        pair<int, int> key {c_i, d_i};
        m_prev_pos[key] = domain->m_pos;
        m_prev_ore[key] = domain->m_ore;
        m_deleted_domains[c_i].push_back(d_i);
        m_delta_e += m_origami_system.unassign_domain(*domain);
    }

    // Delete chain and test acceptance
    m_origami_system.delete_chain(c_i);
    m_deleted_chains.push_back({c_i, c_i_ident});

    return staple_deletion_accepted(c_i_ident);
}

void MetStapleExchangeMCMovetype::select_and_set_growth_point(Domain* growth_domain_new) {

    // Select a different chain, excluding newest addition which is last index
    int c_i_index {m_random_gens.uniform_int(0, m_origami_system.num_staples() - 1)};
    vector<Domain*> selected_chain {m_origami_system.m_domains[c_i_index]};

    // Select a domain of the chain
    int d_i  {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_old {selected_chain[d_i]};

    m_delta_e += set_growth_point(*growth_domain_new, *growth_domain_old);
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
                try {
                    double delta_e {m_origami_system.check_domain_constraints(
                            domain, p_new, o_new)};
                    configs.push_back({p_new, o_new});
                    bfactor *= exp(-delta_e);
                    bfactors.push_back(bfactor);
                }
                catch (ConstraintViolation) {
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
    vector<double> weights {calc_bias(bfactors, domain, configs, p_prev)};

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
    VectorThree o_old {-growth_domain_old.m_ore};
    double delta_e {0};
    try {
        delta_e += m_origami_system.set_domain_config(growth_domain_new,
                growth_domain_old.m_pos, o_old);
        pair<int, int> key {growth_domain_new.m_c, growth_domain_new.m_d};
        m_assigned_domains.push_back(key);
    }
    catch (ConstraintViolation) {
        throw MoveRejection {};
    }

    return delta_e;
}

bool CBMCMovetype::test_cb_acceptance(double new_bias) {
    double ratio {new_bias / m_bias};
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

bool CBStapleExchangeMCMovetype::staple_insertion_accepted(int c_i_ident) {
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    m_bias /= pow(6, staple_length);
    int Ni_new {m_origami_system.num_staples()};

    // Correct for extra states from additional staple domains
    int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
    double extra_states {pow(6, extra_df)};
    double ratio {extra_states / static_cast<double>(Ni_new) * m_bias};

    // Correct for insertion into subset of volume
    m_modifier *= m_insertion_sites / m_origami_system.m_volume;

    // Correct for considering only 1 of staple length ways insertion could occur
    m_modifier *= staple_length;

    return test_acceptance(ratio);
}

bool CBStapleExchangeMCMovetype::staple_deletion_accepted(int c_i_ident) {
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    m_bias /= pow(6, staple_length);
    int Ni_new {m_origami_system.num_staples()};

    // Correct for extra states from additional staple domains
    int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
    double extra_states {pow(6, extra_df)};

    double ratio {(static_cast<double>(Ni_new) + 1) / extra_states / m_bias};

    return test_acceptance(ratio);
}

vector<double> CBStapleExchangeMCMovetype::calc_bias(vector<double> weights,
        Domain*, vector<pair<VectorThree, VectorThree>>& configs, VectorThree) {
    // Calculate rosenbluth weight

    // Set weights of positions involving binding to 0
    for (size_t i {0}; i != configs.size(); i++) {
        if (m_origami_system.m_pos_to_unbound_d.count(configs[i].first) > 0 ) {
            weights[i] = 0;
        }
    }
    double rosenbluth_i {0};
    for (auto weight: weights) {
        rosenbluth_i += weight;
    }
    if (rosenbluth_i == 0) {
        
        // Deadend
        throw MoveRejection {};
    }
    for (size_t i {0}; i != weights.size(); i++) {
        weights[i] /= rosenbluth_i;
    }
    return weights;
}

bool CBStapleExchangeMCMovetype::insert_staple() {

    // Select and add chain of random identity
    int c_i_ident {select_random_staple_identity()};
    int c_i {m_origami_system.add_chain(c_i_ident)};
    m_added_chains.push_back(c_i);

    // Assume that add_chain always adds to end of m_domains
    vector<Domain*> selected_chain {m_origami_system.m_domains.back()};

    // Select growth point on added chain
    int d_i {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_new {selected_chain[d_i]};

    // Select and set growth point on current system
    select_and_set_growth_point(growth_domain_new);

    m_regrow_old = false;
    grow_staple(d_i, selected_chain);

    return staple_insertion_accepted(c_i_ident);
}

bool CBStapleExchangeMCMovetype::delete_staple() {

    // Select random identity and staple of that type
    int c_i_ident {select_random_staple_identity()};
    int c_i {select_random_staple_of_identity(c_i_ident)};

    // The index into m_domains
    int c_i_index {index(m_origami_system.m_chain_indices, c_i)};

    // Find growth piont
    // Also check if staple bound to more than one domain, reject if so
    Domain* growth_domain_new;
    Domain* growth_domain_old;
    vector<Domain*> staple {m_origami_system.m_domains[c_i_index]};
    int bound_domains {0};
    for (auto domain: staple) {
        if (domain->m_state != Occupancy::unbound) {
            bound_domains += 1;
            growth_domain_new = domain;
            growth_domain_old = domain->m_bound_domain;
        }
        if (bound_domains > 1) {
            return false;
        }
    }
    
    // Add stuff to reversions lists and unassign domsins
    int chain_length {static_cast<int>(m_origami_system.m_domains[c_i_index].size())};
    for (int d_i {0}; d_i != chain_length; d_i++) {
        Domain* domain {staple[d_i]};
        pair<int, int> key {c_i, d_i};
        m_old_pos[key] = domain->m_pos;
        m_old_ore[key] = domain->m_ore;
        m_origami_system.unassign_domain(*domain);
    }

    // Regrow chain to calculate Rosenbluth weight
    m_regrow_old = true;
    double delta_e {set_old_growth_point(*growth_domain_new, *growth_domain_old)};
    m_bias *= 6 * exp(-delta_e);
    grow_staple(growth_domain_new->m_d, staple);

    // Prepare reversion list, unassign domains, and delete chain
    m_assigned_domains.clear();
    m_prev_pos = m_old_pos;
    m_prev_ore = m_old_ore;
    for (int d_i {0}; d_i != chain_length; d_i++) {
        Domain* domain {staple[d_i]};
        m_deleted_domains[c_i].push_back(d_i);
        m_origami_system.unassign_domain(*domain);
    }
    m_origami_system.delete_chain(c_i);
    m_deleted_chains.push_back({c_i, c_i_ident});

    return staple_deletion_accepted(c_i_ident);
}

void CBStapleExchangeMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        select_and_set_config(domains, i);
    }
}

void CBStapleExchangeMCMovetype::select_and_set_growth_point(Domain* growth_domain_new) {

    // Select a different chain, excluding newest addition which is last index
    int c_i_index {m_random_gens.uniform_int(0, m_origami_system.num_staples() - 1)};
    vector<Domain*> selected_chain {m_origami_system.m_domains[c_i_index]};

    // Select a domain of the chain
    int d_i  {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_old {selected_chain[d_i]};

    double delta_e {set_growth_point(*growth_domain_new, *growth_domain_old)};
    m_bias *= 6 * exp(-delta_e);
}

bool CBStapleRegrowthMCMovetype::attempt_move() {

    // No staples to regrow
    if (m_origami_system.num_staples() == 0) {

        throw MoveRejection {};
    }

    // Select a staple to regrow
    int c_i_index {m_random_gens.uniform_int(1, m_origami_system.num_staples())};
    vector<Domain*> selected_chain {m_origami_system.m_domains[c_i_index]};

    // Reject if staple is connector
    if (staple_is_connector(selected_chain)) {
        return false;
    }

    // Select growth points on chains
    int growth_d_new {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_new {selected_chain[growth_d_new]};
    Domain* growth_domain_old {select_random_domain()};
    while (growth_domain_old->m_c == growth_domain_new->m_c) {
        growth_domain_old = select_random_domain();
    }

    // Find all bound domains and unassign
    vector<pair<Domain*, Domain*>> bound_domains {};
    for (auto domain: selected_chain) {
        if (domain->m_bound_domain != nullptr) {
            bound_domains.push_back({domain, domain->m_bound_domain});
        }

        pair<int, int> key {domain->m_c, domain->m_d};
        m_prev_pos[key] = domain->m_pos;
        m_prev_ore[key] = domain->m_ore;
        m_modified_domains.push_back(key);

        // Would be double counting to include delta_e in bias
        m_origami_system.unassign_domain(*domain);
    }

    if (bound_domains.empty()) {
        throw OrigamiMisuse {};
    }

    // Grow staple
    m_regrow_old = false;
    double delta_e {0};
    delta_e = set_growth_point(*growth_domain_new, *growth_domain_old);
    m_bias *= exp(-delta_e);
    grow_staple(growth_d_new, selected_chain);

    // Regrow staple in old conformation
    // Save relevant state variables for acceptance testing and resetting
    double new_bias {m_bias};
    m_bias = 1;
    m_modified_domains.clear();
    m_assigned_domains.clear();
    m_old_pos = m_prev_pos;
    m_old_ore = m_prev_ore;

    // Select growth point from previously bound domains
    int bound_domain_index {m_random_gens.uniform_int(0, bound_domains.size() - 1)};
    growth_domain_new = bound_domains[bound_domain_index].first;
    growth_domain_old = bound_domains[bound_domain_index].second;
    
    // Unassign and add to reversion list
    for (auto domain: selected_chain) {
        pair<int, int> key {domain->m_c, domain->m_d};
        m_prev_pos[key] = domain->m_pos;
        m_prev_ore[key] = domain->m_ore;
        m_modified_domains.push_back(key);

        // Would be double counting to include delta_e in bias
        m_origami_system.unassign_domain(*domain);
    }

    // Grow staple
    m_regrow_old = true;
    delta_e = set_old_growth_point(*growth_domain_new, *growth_domain_old);
    m_bias *= exp(-delta_e);
    grow_staple(growth_domain_new->m_d, selected_chain);

    // Test acceptance
    return test_cb_acceptance(new_bias);
}

void CBStapleRegrowthMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        select_and_set_config(domains, i);
    }
}

vector<double> CBStapleRegrowthMCMovetype::calc_bias(vector<double> bfactors,
        Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree) {
    // Calculate rosenbluth weight
    double rosenbluth_i {0};
    for (auto bfactor: bfactors) {
        rosenbluth_i += bfactor;
    }
    if (rosenbluth_i == 0) {
        
        // Deadend
        throw MoveRejection {};
    }
    vector<double> weights {};
    for (auto bfactor: bfactors) {
        weights.push_back(bfactor / rosenbluth_i);
    }
    return weights;
}

bool CTCBScaffoldRegrowthMCMovetype::attempt_move() {

    m_regrow_old = false;
    vector<Domain*> scaffold_domains {select_scaffold_indices()};

    set<int> staples {find_staples_growthpoints_endpoints(scaffold_domains)};
    unordered_map<int, vector<pair<int, VectorThree>>> initial_endpoints {m_active_endpoints};
    
    // Unassign staples except those directly and indirectly bound to external scaffold domains
    vector<Domain*> scaffold_domains_to_unassign(scaffold_domains.begin() + 1,
            scaffold_domains.end());
    unassign_domains(scaffold_domains_to_unassign);
    for (auto c_i: staples) {
        int c_i_index {index(m_origami_system.m_chain_indices, c_i)};
        unassign_domains(m_origami_system.m_domains[c_i_index]);
    }

    // Grow scaffold and staples
    bool growthpoint {m_growthpoints.count(scaffold_domains[0]) > 0};
    if (growthpoint) {
        grow_staple_and_update_endpoints(scaffold_domains, 0);
    }
    grow_chain(scaffold_domains);

    // Regrow in old conformation
    m_regrow_old = true;
    double new_bias {m_bias};
    m_bias = 1;
    m_active_endpoints = initial_endpoints;
    m_modified_domains.clear();
    m_assigned_domains.clear();
    m_old_pos = m_prev_pos;
    m_old_ore = m_prev_ore;

    // Unassign staples except those directly and indirectly bound to external scaffold domains
    unassign_domains(scaffold_domains_to_unassign);
    for (auto c_i: staples) {
        int c_i_index {index(m_origami_system.m_chain_indices, c_i)};
        unassign_domains(m_origami_system.m_domains[c_i_index]);
    }

    // Grow scaffold and staples
    growthpoint = m_growthpoints.count(scaffold_domains[0]) > 0;
    if (growthpoint) {
        grow_staple_and_update_endpoints(scaffold_domains, 0);
    }
    grow_chain(scaffold_domains);

    return test_cb_acceptance(new_bias);
}

vector<Domain*> CTCBScaffoldRegrowthMCMovetype::select_scaffold_indices() {

    // Randomly select endpoints
    vector<Domain*> scaffold {m_origami_system.m_domains[m_origami_system.c_scaffold]};
    Domain* start_domain {scaffold[m_random_gens.uniform_int(0, scaffold.size() - 1)]};
    Domain* end_domain {scaffold[m_random_gens.uniform_int(0, scaffold.size() - 1)]};
    while (start_domain->m_d == end_domain->m_d) {
        end_domain = scaffold[m_random_gens.uniform_int(0, scaffold.size() - 1)];
    }
    
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
        pair<int, VectorThree> endpoint {end_domain->m_d, end_domain->m_pos};
        m_active_endpoints[m_origami_system.c_scaffold].push_back(endpoint);
    }

    // Linear domains
    else {

        // Ensure start domain is 5'
        if (start_domain->m_d > end_domain->m_d) {
            Domain* temp {start_domain};
            start_domain = end_domain;
            end_domain = temp;
        }
        for (int d_i {start_domain->m_d}; d_i != end_domain->m_d + 1; d_i ++) {
            Domain* cur_domain {scaffold[d_i]};
            domains.push_back(cur_domain);
        }

        // If end domain is end of chain, no endpoint
        if (end_domain->m_forward_domain != nullptr) {
            pair<int, VectorThree> endpoint {end_domain->m_d, end_domain->m_pos};
            m_active_endpoints[m_origami_system.c_scaffold].push_back(endpoint);
        }
        else {
        }
    }

    return domains;
}

set<int> CTCBScaffoldRegrowthMCMovetype::find_staples_growthpoints_endpoints(
        vector<Domain*> scaffold_domains) {
    
    unordered_map<Domain*, Domain*> initial_growthpoints = m_growthpoints;
    // Chain ids of staples that have been checked
    set<int> checked_staples {};

    // Staples that are to be regrown
    set<int> regrowth_staples {};
    for (size_t d_i_index {0}; d_i_index != scaffold_domains.size(); d_i_index++) {
        Domain* domain {scaffold_domains[d_i_index]};
        set<int> participating_chains {m_origami_system.c_scaffold};
        vector<pair<Domain*, Domain*>> potential_growthpoints {};
        bool externally_bound {false};
        bool bound {domain->m_bound_domain != nullptr};
        if (bound) {

            // Skip if bound to self
            if (domain->m_bound_domain->m_c == domain->m_c) {
                continue;
            }

            // Skip if already checked
            Domain* bound_domain {domain->m_bound_domain};
            if (find(checked_staples.begin(), checked_staples.end(),
                bound_domain->m_c) != checked_staples.end()) {
                continue;
            }

            potential_growthpoints.push_back({domain, bound_domain});
            scan_staple_topology(bound_domain, participating_chains,
                    potential_growthpoints, scaffold_domains, externally_bound);

            if (externally_bound) {

                // Convert inactive endpoints on scaffold to active
                for (auto c_i: participating_chains) {
                    for (auto endpoint_domain: m_inactive_endpoints[c_i]) {
                        if (endpoint_domain->m_c == m_origami_system.c_scaffold) {
                            pair<int, VectorThree> endpoint {
                                    endpoint_domain->m_d, endpoint_domain->m_pos};

                            m_active_endpoints[m_origami_system.c_scaffold].
                                    push_back(endpoint);
                        }
                    }
                }
            }
            else {
                for (auto growthpoint: potential_growthpoints) {
                    m_growthpoints[growthpoint.first] = growthpoint.second;
                }
                for (auto c_i: participating_chains) {
                    if (c_i == 0) {
                        continue;
                    }
                    regrowth_staples.insert(c_i);
                }
            }

            // Add to checked staples
            for (auto c_i: participating_chains) {
                checked_staples.insert(c_i);
            }
        }
    }

    // Remove endpoints on first scaffold domain
    vector<int> endpoints_to_erase {};
    for (size_t j {0}; j != m_active_endpoints[m_origami_system.c_scaffold].size(); j++) {
        pair<int, VectorThree> endpoint {m_active_endpoints[m_origami_system.c_scaffold][j]};
        if (endpoint.first == scaffold_domains[0]->m_d) {
            endpoints_to_erase.push_back(j);
        }
    }
    std::reverse(endpoints_to_erase.begin(), endpoints_to_erase.end());
    for (auto j: endpoints_to_erase) {
        m_active_endpoints[m_origami_system.c_scaffold].erase(
                m_active_endpoints[m_origami_system.c_scaffold].begin() + j);
    }

    return regrowth_staples;
}

void CTCBScaffoldRegrowthMCMovetype::unassign_domains(vector<Domain*> domains) {
    for (auto domain: domains) {
        pair<int, int> key {domain->m_c, domain->m_d};
        m_prev_pos[key] = domain->m_pos;
        m_prev_ore[key] = domain->m_ore;
        m_modified_domains.push_back(key);

        // Would be double counting to include delta_e in bias
        m_origami_system.unassign_domain(*domain);
    }
}

void CTCBScaffoldRegrowthMCMovetype::scan_staple_topology(
        Domain* domain,
        set<int>& participating_chains,
        vector<pair<Domain*, Domain*>>& potential_growthpoints,
        vector<Domain*>& scaffold_domains,
        bool& externally_bound) {

    participating_chains.insert(domain->m_c);
    int c_i_index {index(m_origami_system.m_chain_indices, domain->m_c)};
    vector<Domain*> staple {m_origami_system.m_domains[c_i_index]};
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

            // Check if bound domain already in progress
            if (participating_chains.count(bound_domain->m_c) > 0) {

                // Check if externally bound
                bool domain_on_scaffold {bound_domain->m_c == m_origami_system.c_scaffold};
                bool domain_in_range {find(scaffold_domains.begin(), scaffold_domains.end(),
                        bound_domain) != scaffold_domains.end()};
                if (domain_on_scaffold) {
                    if (domain_in_range) {
                        m_inactive_endpoints[bound_domain->m_c].push_back(bound_domain);
                    }
                    else {
                        externally_bound = true;
                    }
                }
            }
            else {
                potential_growthpoints.push_back({cur_domain, bound_domain});
                scan_staple_topology(bound_domain, participating_chains,
                        potential_growthpoints, scaffold_domains, externally_bound);
            }
        }
    }
}

void CTCBScaffoldRegrowthMCMovetype::grow_chain(vector<Domain*> domains) {
    int c_i {domains[0]->m_c};
    for (size_t i {1}; i != domains.size(); i++) {
        select_and_set_config(domains, i);
    
        // Remove endpoint if reached
        vector<int> endpoints_to_erase {};
        for (size_t j {0}; j != m_active_endpoints[c_i].size(); j++) {
            pair<int, VectorThree> endpoint {m_active_endpoints[c_i][j]};
            if (endpoint.first == domains[i]->m_d) {
                endpoints_to_erase.push_back(j);
            }
        }
        std::reverse(endpoints_to_erase.begin(), endpoints_to_erase.end());
        for (auto j: endpoints_to_erase) {
            m_active_endpoints[c_i].erase(m_active_endpoints[c_i].begin() + j);
        }

        // Grow staple if growth point
        bool growthpoint {m_growthpoints.count(domains[i]) > 0};
        if (growthpoint) {
            grow_staple_and_update_endpoints(domains, i);
        }
    }
}

void CTCBScaffoldRegrowthMCMovetype::grow_staple_and_update_endpoints(
        vector<Domain*> domains, int i) {
    Domain* growth_domain_old {domains[i]};
    Domain* growth_domain_new {m_growthpoints[domains[i]]};
    int c_i {growth_domain_new->m_c};
    double delta_e;
    if (m_regrow_old) {
        delta_e = set_old_growth_point(*growth_domain_new, *growth_domain_old);
    }
    else {
        delta_e = set_growth_point(*growth_domain_new, *growth_domain_old);
    }
    m_bias *= exp(-delta_e);
    int c_i_index {index(m_origami_system.m_chain_indices,
            growth_domain_new->m_c)};
    vector<Domain*> staple {m_origami_system.m_domains[c_i_index]};
    grow_staple(growth_domain_new->m_d, staple);

    // Copy inactive endpoints to active endpoints for finished chain
    bool inactive_endpoints_present {m_inactive_endpoints.count(c_i) > 0};
    if (inactive_endpoints_present) {
        for (auto domain: m_inactive_endpoints[c_i]) {
            pair<int, VectorThree> endpoint {domain->m_d, domain->m_pos};
            m_active_endpoints[domain->m_c].push_back(endpoint);
        }
    }
}

vector<double> CTCBScaffoldRegrowthMCMovetype::calc_bias(
        vector<double> weights,
        Domain* domain,
        vector<pair<VectorThree, VectorThree>>& configs,
        VectorThree p_prev) {

    // Set weights of positions involving non-self binding to 0 unless endpoint reached
    for (size_t i {0}; i != configs.size(); i++) {
        if (m_origami_system.m_pos_to_unbound_d.count(configs[i].first) > 0 ) {
            if (m_origami_system.unbound_domain_at(configs[i].first)->m_d == domain->m_d) {
                continue;
            }
            for (auto endpoint: m_active_endpoints[domain->m_c]) {
                bool endpoint_domain_reached {domain->m_d == endpoint.first};
                bool endpoint_pos_reached {configs[i].first == endpoint.second};
                if (endpoint_domain_reached and endpoint_pos_reached) {
                    continue;
                }
                else {
                    weights[i] = 0;
                }
            }
        }
    }

    // Bias weights with number of walks
    for (size_t i {0}; i != configs.size(); i++) {
        VectorThree start_pos {configs[i].first};
        double prod_num_walks {1};
        for (auto endpoint: m_active_endpoints[domain->m_c]) {
            int steps {endpoint.first - domain->m_d};
            VectorThree endpoint_p {endpoint.second};
            prod_num_walks *= m_ideal_random_walks.num_walks(start_pos,
                    endpoint_p, steps);
        }
        weights[i] *= prod_num_walks;
    }

    // Calculate number of walks for previous position
    double prod_num_walks {1};
    for (auto endpoint: m_active_endpoints[domain->m_c]) {
        int steps {endpoint.first - domain->m_d + 1};
        VectorThree endpoint_p {endpoint.second};
        prod_num_walks *= m_ideal_random_walks.num_walks(p_prev,
                endpoint_p, steps);
    }

    // Modified Rosenbluth
    double weights_sum {0};

    for (auto weight: weights) {
        weights_sum += weight;
    }

    // Check for deadend
    if (weights_sum == 0) {
        throw MoveRejection {};
    }

    double bias {weights_sum / prod_num_walks};
    m_bias *= bias;

    // Normalize
    for (size_t i {0}; i != weights.size(); i++) {
        weights[i] /= weights_sum;
    }

    return weights;
}

template<typename T>
unique_ptr<MCMovetype> Movetypes::movetype_constructor(OrigamiSystem& origami_system, RandomGens& random_gens) {
    return unique_ptr<MCMovetype> {new T {origami_system, random_gens}};
}
