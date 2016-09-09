// movetypes.cpp

#include <algorithm>
#include <random>

#include "utility.h"
#include "random_gens.h"
#include "movetypes.h"

using std::fmin;

using namespace Movetypes;
using namespace Utility;
using namespace RandomGen;

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
            try{
            m_origami_system.set_checked_domain_config(*domain, pos, ore);
    } catch (ConstraintViolation) {cout << "e\n";}
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
        try{
        m_origami_system.set_checked_domain_config(*domain, pos, ore);
    } catch (ConstraintViolation) {cout << "f\n";}
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
    int d_i_index {m_random_gens.uniform_int(0, m_origami_system.m_num_domains - 1)};
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
    int old_num_bound_domains {m_origami_system.m_num_domains};

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
    int overcount {old_num_bound_domains - m_origami_system.m_num_domains};
    m_modifier /= overcount;
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

bool StapleExchangeMCMovetype::attempt_move() {

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

bool StapleExchangeMCMovetype::staple_insertion_accepted(int c_i_ident) {
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

bool StapleExchangeMCMovetype::staple_deletion_accepted(int c_i_ident) {
    double boltz_factor {exp(-m_delta_e)};
    int Ni_new {m_origami_system.num_staples()};

    // Correct for extra states from additional staple domains
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
    int extra_states {static_cast<int>(pow(6, extra_df))};

    double ratio {(Ni_new + 1) / extra_states * boltz_factor};

    return test_acceptance(ratio);
}

bool StapleExchangeMCMovetype::insert_staple() {

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

bool StapleExchangeMCMovetype::delete_staple() {

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

void StapleExchangeMCMovetype::select_and_set_growth_point(Domain* growth_domain_new) {

    // Select a different chain, excluding newest addition which is last index
    int c_i_index {m_random_gens.uniform_int(0, m_origami_system.num_staples() - 1)};
    vector<Domain*> selected_chain {m_origami_system.m_domains[c_i_index]};

    // Select a domain of the chain
    int d_i  {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_old {selected_chain[d_i]};

    m_delta_e += set_growth_point(*growth_domain_new, *growth_domain_old);
}

bool CBStapleRegrowthMCMovetype::attempt_move() {

    // No staples to regrow
    if (m_origami_system.num_staples() == 0) {

        throw MoveRejection {};
    }

    // Select a staple to regrow and a growth point on it
    int c_i_index {m_random_gens.uniform_int(1, m_origami_system.num_staples())};
    vector<Domain*> selected_chain {m_origami_system.m_domains[c_i_index]};
    int growth_d_new {m_random_gens.uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_new {selected_chain[growth_d_new]};

    // Select a growth point on one of the other chains
    // THIS IS WRONG
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
    delta_e = set_growth_point(*growth_domain_new, *growth_domain_old);
    m_bias *= exp(-delta_e);
    grow_staple(growth_domain_new->m_d, selected_chain);

    // Test acceptance
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

    return false;
}

void CBStapleRegrowthMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        Domain* domain {domains[i]};
        Domain* prev_domain {domains[i - 1]};
        VectorThree p_prev {prev_domain->m_pos};

        // Calculate weights
        vector<pair<VectorThree, VectorThree>> configs {};
        vector<double> bfactors {};
        calculate_biases(*domain, p_prev, configs, bfactors);
        vector<double> weights {calc_bias(bfactors)};

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
}

void CBStapleRegrowthMCMovetype::calculate_biases(
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

vector<double> CBStapleRegrowthMCMovetype::calc_bias(vector<double> bfactors) {
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

void CBStapleRegrowthMCMovetype::select_and_set_new_config(
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
        }
    }

    // Set config on domain
    VectorThree p_new {config.first};
    VectorThree o_new {config.second};
    if (o_new == VectorThree {0, 0, 0}) {
        o_new = select_random_orientation();
        try{
        m_origami_system.set_checked_domain_config(domain, p_new, o_new);
    } catch (ConstraintViolation) {cout << "g\n";}
    }
    else {
        try{
        m_origami_system.set_checked_domain_config(domain, p_new, o_new);
    } catch (ConstraintViolation) {cout << "h\n";}
    }
}

void CBStapleRegrowthMCMovetype::select_and_set_old_config(Domain& domain) {
    pair<int, int> key {domain.m_c, domain.m_d};
    VectorThree p_old {m_old_pos[key]};
    VectorThree o_old {m_old_ore[key]};
    try{
    m_origami_system.set_checked_domain_config(domain, p_old, o_old);
    } catch (ConstraintViolation) {cout << "d\n";}
}

template<typename T>
unique_ptr<MCMovetype> Movetypes::movetype_constructor(OrigamiSystem& origami_system, RandomGens& random_gens) {
    return unique_ptr<MCMovetype> {new T {origami_system, random_gens}};
}

