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
    std::uniform_int_distribution<int> dist(0, m_origami_system.m_domains.size());
    auto gen_chain = bind(dist, random_engine);
    int c_i_index {gen_chain()};
    std::uniform_int_distribution<int> dist_domain(0, m_origami_system.m_domains[c_i_index].size());
    auto gen_domain = bind(dist_domain, random_engine);
    Domain* domain {m_origami_system.m_domains[c_i_index][gen_domain()]};

    return domain;
}

int MCMovetype::select_random_staple_identity() {

    // Randomly select staple identity and then staple of that identity
    std::uniform_int_distribution<int> dist(1, m_origami_system.m_identities.size());
    auto gen = bind(dist, random_engine);
    int c_i_ident {gen()};

    return c_i_ident;
}

int MCMovetype::select_random_staple_of_identity(int c_i_ident) {
    size_t staples {m_origami_system.m_identity_to_index[c_i_ident].size()};
    if (staples == 0) {
        throw MoveRejection {};
    }
    std::uniform_int_distribution<int> dist2(0, staples);
    auto gen2 = bind(dist2, random_engine);
    int c_i {m_origami_system.m_identity_to_index[c_i_ident][gen2()]};

    return c_i;
}

double MCMovetype::set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old) {
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

void MCMovetype::calc_overcount(Domain& domain) {
    if (domain.m_bound_domain != nullptr) {
        size_t chain_length {m_origami_system.m_identities[domain.m_c_ident].size()};
    }
}

template<typename T>
unique_ptr<MCMovetype> Movetypes::movetype_constructor(OrigamiSystem origami_system) {
    return new T {origami_system};
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
    VectorThree vec {vectors[gen_uniform_vector_i()]};
    return p_prev + vec;
}

VectorThree Movetypes::select_random_orientation() {
    VectorThree vec {vectors[gen_uniform_vector_i()]};
    return vec;
}

bool OrientationRotationMCMovetype::attempt_move() {

    bool accept {false};
    
    // Select random chain and domain
    Domain* domain {select_random_domain()};

    // Reject if in bound state
    if (domain->m_state == Occupancy::bound) {
        accept = false;
        return accept;
    }
    else {
        accept = true;
    }

    // Select random orientation and update
    VectorThree o_new {select_random_orientation()};
    m_origami_system.set_domain_orientation(*domain, o_new);

    return accept;
}

void MetMCMovetype::grow_staple(int d_i_index, vector<Domain*> selected_chain) {
    auto first_iter {selected_chain.begin() + d_i_index};
    auto last_iter {selected_chain.end() + selected_chain.size()};
    vector<Domain*> domains_three_prime {first_iter, last_iter};
    grow_chain(domains_three_prime);
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
        calc_overcount(*domain);
    }
}

bool StapleExchangeMCMovetype::attempt_move() {
    bool accept {false};
    if (gen_uniform_real() < 0.5) {
        accept = insert_staple();
    }
    else {
        accept = delete_staple();
    }
    
    return accept;
}

bool StapleExchangeMCMovetype::staple_insertion_accepted(int c_i_ident) {
    double boltz_factor {exp(-m_delta_e / m_origami_system.m_temp)};
    int Ni_new {m_origami_system.num_staples()};

    // Correct for extra states from additional staple domains
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    int extra_states {static_cast<int>(
            pow(6, 2 * staple_length - 1 - preconstrained_df))};
    double ratio {extra_states / Ni_new * boltz_factor};

    // Correct for insertion into subset of volume
    m_modifier *= m_insertion_sites / m_origami_system.m_volume;

    // Correct for considering only 1 of staple length ways insertion could occur
    m_modifier *= staple_length;

    return test_acceptance(ratio, m_modifier);
}

bool StapleExchangeMCMovetype::staple_deletion_accepted(int c_i_ident) {
    double boltz_factor {exp(-m_delta_e / m_origami_system.m_temp)};
    int Ni_new {m_origami_system.num_staples()};

    // Correct for extra states from additional staple domains
    size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
    int extra_states {static_cast<int>(
            pow(6, 2 * staple_length - 1 - preconstrained_df))};
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

    // Select growth points and attempt binding
    std::uniform_int_distribution<int> dist2(1, selected_chain.size());
    auto gen2 = bind(dist2, random_engine);
    int d_i_index {gen2()};
    Domain* growth_domain_new {selected_chain[d_i_index]};
    select_and_set_growth_point(growth_domain_new);

    grow_staple(d_i_index, selected_chain);

    return staple_insertion_accepted(c_i_ident);
}

bool StapleExchangeMCMovetype::delete_staple() {

    int c_i_ident {select_random_staple_identity()};
    int c_i {select_random_staple_of_identity(c_i_ident)};
    int c_i_index {index(m_origami_system.m_chain_indices, c_i)};
    
    // Add stuff to reversions lists and unassign domsins
    for (auto domain: m_origami_system.m_domains[c_i_index]) {
        m_prev_pos[domain] = domain->m_pos;
        m_prev_ore[domain] = domain->m_ore;
        m_deleted_domains.push_back(domain);
        m_delta_e += m_origami_system.unassign_domain(*domain);
    }
    m_origami_system.delete_chain(c_i);
    m_deleted_chains.push_back({c_i, c_i_ident});

    return staple_deletion_accepted(c_i_ident);
}

void StapleExchangeMCMovetype::select_and_set_growth_point(Domain* growth_domain_new) {
    std::uniform_int_distribution<int> dist(1, m_origami_system.num_staples());
    auto gen = bind(dist, random_engine);
    int c_i_index {gen()};
    vector<Domain*> selected_chain {m_origami_system.m_domains[c_i_index]};
    std::uniform_int_distribution<int> dist2(1, m_origami_system.num_staples());
    auto gen2 = bind(dist2, random_engine);
    Domain* growth_domain_old {selected_chain[gen2()]};
    m_delta_e += set_growth_point(*growth_domain_new, *growth_domain_old);
}

bool CBStapleRegrowthMCMovetype::attempt_move() {

    // Pick section of scaffold to regrow and determine
    
}
