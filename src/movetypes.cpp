// movetypes.cpp

#include <algorithm>
#include <random>
#include <set>
#include <cmath>

#include "utility.h"
#include "movetypes.h"
#include "ideal_random_walk.h"

using std::fmin;
using std::set;
using std::find;
using std::pow;

using namespace Movetypes;
using namespace Utility;

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

bool MCMovetype::test_acceptance(long double p_ratio) {
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
    //// Set growth point with new orientation
    //VectorThree o_new {select_random_orientation()};
    // Set growth point with complementary orientation
    double delta_e {0};
    VectorThree o_new {-growth_domain_old.m_ore};
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
    int old_num_bound_domains {m_origami_system.num_bound_domain_pairs() -
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
    int overcount {m_origami_system.num_bound_domain_pairs() -
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

template<typename T>
unique_ptr<MCMovetype> Movetypes::movetype_constructor(
        OrigamiSystem& origami_system,
        RandomGens& random_gens,
        IdealRandomWalks& ideal_random_walks,
        InputParameters& params) {
    return unique_ptr<MCMovetype> {new T {origami_system, random_gens,
            ideal_random_walks, params}};
}
