// top_constraint_points.cpp

#include <set>
#include <cmath>

#include "utility.h"
#include "top_constraint_points.h"

using std::fmin;
using std::set;
using std::find;

namespace TopConstraintPoints {

    using namespace Utility;

    Constraintpoints::Constraintpoints(
            OrigamiSystem& origami_system,
            IdealRandomWalks& ideal_random_walks) :
            m_origami_system {origami_system},
            m_ideal_random_walks {ideal_random_walks} {
    }

    void Constraintpoints::reset_internal() {
        m_scaffold_domains.clear();
        m_regrowth_staples.clear();
        m_growthpoints.clear();
        m_active_endpoints.clear();
        m_initial_active_endpoints.clear();
        m_inactive_endpoints.clear();
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

    long double Constraintpoints::calc_num_walks_prod(
            Domain* domain,
            VectorThree pos,
            vector<Domain*> domains,
            // Direction of growth (1 for forward, -1 for backward)
            int dir,
            // For calculating previous position endpoints
            int step_offset) {

        long double prod_num_walks {1};
        for (auto endpoint: m_active_endpoints[domain->m_c]) {

            // Check if endpoint in current stretch of domains being regrown (
            // staples are grown out from their growth point in two seperate
            // directions)
            int endpoint_d_i {endpoint.first};
            int first_d_i {domains.front()->m_d};
            int last_d_i {domains.back()->m_d};
            bool apply_constraint {false};
            // Why is it not that condition?
            if (domain->m_c != m_origami_system.c_scaffold) {
                if ((endpoint_d_i >= first_d_i and endpoint_d_i <= last_d_i) or
                        (endpoint_d_i >= last_d_i and endpoint_d_i <= first_d_i)) {
                    apply_constraint = true;
                }
            }
            else {
                if ((endpoint_d_i >= first_d_i and endpoint_d_i <= (last_d_i + 1)) or
                        (endpoint_d_i >= (last_d_i - 1) and endpoint_d_i <= (first_d_i))) {
                    apply_constraint = true;
                }
            }
            if (apply_constraint) {
                int steps {calc_remaining_steps(endpoint_d_i, domain, dir,
                        step_offset)};
                VectorThree endpoint_p {endpoint.second};
                prod_num_walks *= m_ideal_random_walks.num_walks(pos, endpoint_p,
                        steps);
            }
        }

        return prod_num_walks;
    }

    int Constraintpoints::calc_remaining_steps(int endpoint_d_i, Domain* domain,
            int dir, int step_offset) {
        int steps;
        if (m_origami_system.m_cyclic and domain->m_c == m_origami_system.c_scaffold) {
            if (dir > 0 and endpoint_d_i < domain->m_d) {
                steps = m_origami_system.get_chain(0).size() + endpoint_d_i - domain->m_d;
            }
            else if (dir < 0 and endpoint_d_i > domain->m_d) {
                steps = domain->m_d + m_origami_system.get_chain(0).size() - endpoint_d_i;
            }
            else if (dir == 0 or endpoint_d_i == domain->m_d) {
                steps = 0;
            }
            else {
                steps = abs(endpoint_d_i - domain->m_d);
            }
        }
        else {
            steps = abs(endpoint_d_i - domain->m_d);
        }
        steps += step_offset;

        return steps;
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

        // Remove endpoints on first scaffold domain unless cyclic and regrowing whole thing
        if (not m_origami_system.m_cyclic and scaffold_domains.size() !=
                m_origami_system.get_chain(0).size()) {
            remove_active_endpoint(scaffold_domains[0]);
        }
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
}
