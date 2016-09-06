// movetypes.cpp

#include <algorithm>
#include <random>

#include "utility.h"
#include "movetypes.h"

using std::fmin;

using namespace Movetypes;
using namespace Utility;

void MCMovetype::reset_origami() {
    
    // Have to deal with added domains, modified domains, and deleted domains

    // Unassign assigned domains
    for (auto domain: m_assigned_domains) {
        m_origami_system.unassign_domain(*domain);
    }

    // Delete chains that were added
    for (auto c_i: m_added_chains) {
        m_origami_system.delete_chain(c_i);
    }

    // Re-add chains that were deleted revert domains to previous positions
    for (auto index_identity: m_deleted_chains) {
        m_origami_system.add_chain(index_identity.first, index_identity.second);
        for (auto domain: m_deleted_domains) {
            VectorThree pos {m_prev_pos[domain]};
            VectorThree ore {m_prev_ore[domain]};
            m_origami_system.set_checked_domain_config(*domain, pos, ore);
        }
    }

    // Revert modified domains to previous positions
    for (auto domain: m_modified_domains) {
        VectorThree pos {m_prev_pos[domain]};
        VectorThree ore {m_prev_ore[domain]};
        m_origami_system.set_checked_domain_config(*domain, pos, ore);
    }

}

Domain* MCMovetype::select_random_domain() {
    int c_i_index {gen_uniform_int(0, m_origami_system.m_domains.size() - 1)};
    int d_i_index {gen_uniform_int(0, m_origami_system.m_domains[c_i_index].size() - 1)};
    Domain* domain {m_origami_system.m_domains[c_i_index][d_i_index]};

    return domain;
}

int MCMovetype::select_random_staple_identity() {
    return gen_uniform_int(1, m_origami_system.m_identities.size() - 1);
}

int MCMovetype::select_random_staple_of_identity(int c_i_ident) {
    size_t staples {m_origami_system.m_identity_to_index[c_i_ident].size()};
    if (staples == 0) {
        throw MoveRejection {};
    }
    int staple_i {gen_uniform_int(0, staples - 1)};
    int c_i {m_origami_system.m_identity_to_index[c_i_ident][staple_i]};

    return c_i;
}

double MCMovetype::set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old) {
    // Set growth point with new orientation
    VectorThree o_new = select_random_orientation();
    double delta_e {0};
    try {
        delta_e += m_origami_system.set_domain_config(growth_domain_new,
                growth_domain_old.m_pos, o_new);
        m_assigned_domains.push_back(&growth_domain_new);
    }
    catch (ConstraintViolation) {
        throw MoveRejection {};
    }

    return delta_e;
}

template<typename T>
unique_ptr<MCMovetype> Movetypes::movetype_constructor(OrigamiSystem& origami_system) {
    return unique_ptr<MCMovetype> {new T {origami_system}};
}

bool Movetypes::test_acceptance(double p_ratio, double modifier) {
    double p_accept = fmin(1, p_ratio) * modifier;
    bool accept;
    if (p_accept == 1) {
        accept = true;
    }
    else {
        if (p_accept > gen_uniform_real()) {
            accept = true;
        }
        else {
            accept = false;
        }
    }

    return accept;
}

VectorThree Movetypes::select_random_position(VectorThree p_prev) {
    VectorThree vec {vectors[gen_uniform_int(0, 5)]};
    return p_prev + vec;
}

VectorThree Movetypes::select_random_orientation() {
    VectorThree vec {vectors[gen_uniform_int(0, 5)]};
    return vec;
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

void MetMCMovetype::grow_staple(int d_i_index, vector<Domain*> selected_chain) {
    int old_num_bound_domains {m_origami_system.m_num_domains};

    // Grow in three prime direction
    auto first_iter3 {selected_chain.begin() + d_i_index};
    auto last_iter3 {selected_chain.end() + selected_chain.size()};
    vector<Domain*> domains_three_prime {first_iter3, last_iter3};
    grow_chain(domains_three_prime);

    // Grow in five prime direction
    auto first_iter5 {selected_chain.begin()};
    auto last_iter5 {selected_chain.begin() + d_i_index};
    vector<Domain*> domains_five_prime {first_iter5, last_iter5};
    std::reverse(domains_five_prime.begin(), domains_five_prime.end());
    grow_chain(domains_three_prime);

    // Overcount correction
    int overcount {old_num_bound_domains - m_origami_system.m_num_domains};
    m_modifier /= overcount;
}

void MetMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        Domain* domain {domains[i]};
        Domain* prev_domain {domains[i - 1]};
        VectorThree new_p {select_random_position(prev_domain->m_pos)};
        VectorThree new_o {select_random_orientation()};
        try {
            m_origami_system.set_domain_config(*domain, new_p, new_o);
            m_assigned_domains.push_back(domain);
        }
        catch (ConstraintViolation) {
            throw MoveRejection {};
        }
    }
}

bool StapleExchangeMCMovetype::attempt_move() {

    // Select inertion or deletion with equal frequency
    bool accept;
    if (gen_uniform_real() < 0.5) {
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

    return test_acceptance(ratio, m_modifier);
}

bool StapleExchangeMCMovetype::staple_deletion_accepted(int c_i_ident) {
    double boltz_factor {exp(-m_delta_e)};
    int Ni_new {m_origami_system.num_staples()};

    // Correct for extra states from additional staple domains
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
    int extra_states {static_cast<int>(pow(6, extra_df))};

    double ratio {(Ni_new + 1) / extra_states * boltz_factor};

    return test_acceptance(ratio, m_modifier);
}

bool StapleExchangeMCMovetype::insert_staple() {
    
    // Select and add chain of random identity
    int c_i_ident {select_random_staple_identity()};
    int c_i {m_origami_system.add_chain(c_i_ident)};
    m_added_chains.push_back(c_i);

    // Assume that add_chain always adds to end of m_domains
    vector<Domain*> selected_chain {m_origami_system.m_domains.back()};

    // Add domains of given chain to reversion tracker list
    for (auto domain: selected_chain) {
        m_added_domains.push_back(domain);
    }

    // Select growth point on added chain
    int d_i_index {gen_uniform_int(1, selected_chain.size() - 1)};
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
    for (auto domain: m_origami_system.m_domains[c_i_index]) {
        m_prev_pos[domain] = domain->m_pos;
        m_prev_ore[domain] = domain->m_ore;
        m_deleted_domains.push_back(domain);
        m_delta_e += m_origami_system.unassign_domain(*domain);
    }

    // Delete chain and test acceptance
    m_origami_system.delete_chain(c_i);
    m_deleted_chains.push_back({c_i, c_i_ident});

    return staple_deletion_accepted(c_i_ident);
}

void StapleExchangeMCMovetype::select_and_set_growth_point(Domain* growth_domain_new) {

    // Select a different chain, excluding newest addition which is last index
    int c_i_index {gen_uniform_int(0, m_origami_system.num_staples() - 1)};
    vector<Domain*> selected_chain {m_origami_system.m_domains[c_i_index]};

    // Select a domain of the chain
    int d_i_index  {gen_uniform_int(0, selected_chain.size() - 1)};
    Domain* growth_domain_old {selected_chain[d_i_index]};

    m_delta_e += set_growth_point(*growth_domain_new, *growth_domain_old);
}

bool CBStapleRegrowthMCMovetype::attempt_move() {

    // Pick section of scaffold to regrow and determine
    return false;
    
}
