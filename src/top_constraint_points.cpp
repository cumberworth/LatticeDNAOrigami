// top_constraint_points.cpp

#include <cmath>

#include "utility.h"
#include "top_constraint_points.h"

namespace topConstraintPoints {

    using std::fmin;
    using std::set;
    using std::find;

    bool domain_included(vector<Domain*> domains, Domain* d) {
        return find(domains.begin(), domains.end(), d) !=
            domains.end();
    }

    bool chain_included(vector<int> staples, int staple) {
        return find(staples.begin(), staples.end(), staple) !=
        staples.end();
    }

    bool chain_included(set<int> staples, int staple) {
        return staples.find(staple) != staples.end();
    }

    StapleNetwork::StapleNetwork(OrigamiSystem& origami):
        m_origami {origami} {
    }

    void StapleNetwork::scan_network(
            Domain* d,
            vector<int> ex_staples) {
        
        clear_network();
        m_ex_staples = ex_staples;

        // Adding the scaffold chain made checks easier
        m_net_cs.insert(m_origami.c_scaffold);

        // Start recursive scan of network
        scan_staple_topology(d);
    }

    set<int> StapleNetwork::get_participating_chains() {
        return m_net_cs;
    }

    vector<pair<Domain*, Domain*>> StapleNetwork::get_potential_growthpoints() {
        return m_pot_gps;
    }

    vector<pair<Domain*, Domain*>> StapleNetwork::get_potential_inactive_endpoints() {
        return m_pot_iaes;
    }

    unordered_map<Domain*, int> StapleNetwork::get_staple_to_segs_map() {
        return m_segs;
    }

    bool StapleNetwork::externally_bound() {
        return m_external;
    }

    void StapleNetwork::clear_network() {
        m_external = false;
        m_net_cs.clear();
        m_pot_gps.clear();
        m_pot_ds.clear();
            
    }

    void StapleNetwork::scan_staple_topology(Domain* growth_d) {

        int ci {growth_d->m_c}; // Staple chain index
        m_net_cs.insert(ci);
        auto ds {make_staple_stack(growth_d, ci)}; // Domains to check
        for (auto d: ds) {
            m_pot_ds.push_back(d);
            Domain* bd {d->m_bound_domain}; // Domain bound to current
            bool domain_bound {bd != nullptr};

            // Skip if not bound or bound to self
            if (not domain_bound or (bd->m_c == ci)) {
                continue;
            }
            int bd_ci {bd->m_c}; // Bound domain chain index

            /*
             * If the staple has not be added to the network, recurse into the
             * bound staple. Otherwise consider adding potential endpoints
             */
            bool bs_excluded {chain_included(m_ex_staples, bd_ci)};
            bool in_network {chain_included(m_net_cs, bd_ci)};
            if (in_network) {

                /*
                 * If the bound domain is on the scaffold and in the range being
                 * regrown or on another staple, and nothing is excluded the
                 * domans are potential endpoints. Otherwise the network is
                 * externally bound
                 */
                bool d_excluded {false};
                if (bd_ci == 0) {
                    m_external = not domain_included(m_scaffold_ds, bd);
                }
                else {
                    d_excluded = chain_included(m_ex_staples, ci);
                }
                if (not m_external and not (bs_excluded and d_excluded)) {
                    add_potential_inactive_endpoint(d, bd);

                }
            }
            else {

                /*
                 * Add domain to potential growthpoints if bound domain if not
                 * on an excluded staple
                 */
                if (not bs_excluded) {
                    m_pot_gps.push_back({d, bd});
                }
                scan_staple_topology(bd);
            }
        }
    }

    vector<Domain*> StapleNetwork::make_staple_stack(Domain* d, int ci) {
        vector<Domain*> ds {};
        vector<Domain*> staple {m_origami.get_chain(ci)};

        // Iterate through domains in three prime direction
        int seg {0};
        m_segs[d] = seg;
        for (size_t di = d->m_d + 1; di != staple.size(); di++) {
            Domain* d_loop {staple[di]};
            ds.push_back(d_loop);
            m_segs[d_loop] = seg;
        }

        // Add domains in five prime direction
        for (int di {d->m_d - 1}; di != -1; di--) {
            ds.push_back(staple[di]);
        }
        
        return ds;
    }

    void StapleNetwork::add_potential_inactive_endpoint(Domain* d, Domain* bd) {
        if (domain_included(m_pot_ds, bd)) {
            m_pot_iaes.push_back({bd, d});
        }
        else {
            m_pot_iaes.push_back({d, bd});
        }
    }

    Constraintpoints::Constraintpoints(
            OrigamiSystem& origami_system,
            IdealRandomWalks& ideal_random_walks) :
            m_origami_system {origami_system},
            m_ideal_random_walks {ideal_random_walks},
            m_staple_network {origami_system} {
    }

    void Constraintpoints::reset_internal() {
        m_scaffold_domains.clear();
        m_regrowth_staples.clear();
        m_growthpoints.clear();
        m_active_endpoints.clear();
        m_initial_active_endpoints.clear();
        m_inactive_endpoints.clear();
        m_d_stack.clear();
        m_segs.clear();
    }

    void Constraintpoints::calculate_constraintpoints(
            vector<Domain*> scaffold_domains,
            vector<int> excluded_staples) {

        find_growthpoints_endpoints(scaffold_domains, excluded_staples, 0);

        // Save initial active endpionts and remaining steps for regrowing old
        m_initial_active_endpoints = m_active_endpoints;
    }

    void Constraintpoints::calculate_constraintpoints(
            vector<vector<Domain*>> scaffold_segments,
            vector<int> excluded_staples) {

        int seg {0};
        for (auto scaffold_domains: scaffold_segments) {
            find_growthpoints_endpoints(scaffold_domains, excluded_staples, seg);
            seg++;
        }

        // Save initial active endpionts and remaining steps for regrowing old
        m_initial_active_endpoints = m_active_endpoints;
    }

    set<int> Constraintpoints::staples_to_be_regrown() {
        return m_regrowth_staples;
    }

    vector<Domain*> Constraintpoints::domains_to_be_regrown() {
        return m_d_stack;
    }

    bool Constraintpoints::is_growthpoint(Domain* domain) {
        return (m_growthpoints.count(domain) > 0);
    }

    void Constraintpoints::add_active_endpoint(
            Domain* d,
            VectorThree pos) {

        auto seg = m_segs[d];
        pair<int, int> key {d->m_c, seg};
        m_active_endpoints[key].push_back({d->m_d, pos});
    }

    void Constraintpoints::add_active_endpoint(
            Domain* domain,
            VectorThree pos,
            int seg) {

        pair<int, int> key {domain->m_c, seg};
        m_active_endpoints[key].push_back({domain->m_d, pos});
    }

    void Constraintpoints::reset_active_endpoints() {
        m_active_endpoints = m_initial_active_endpoints;
    }

    void Constraintpoints::remove_active_endpoint(Domain* domain) {

        // Remove endpoint if reached
        int c_i {domain->m_c};
        auto seg = m_segs[domain];
        pair<int, int> key {c_i, seg};
        vector<int> endpoints_to_erase {};
        for (size_t j {0}; j != m_active_endpoints[key].size(); j++) {
            pair<int, VectorThree> endpoint {m_active_endpoints[key][j]};
            if (endpoint.first == domain->m_d) {
                endpoints_to_erase.push_back(j);
            }
        }

        // Removing in reverse means the indices are not changing as I do it
        std::reverse(endpoints_to_erase.begin(), endpoints_to_erase.end());
        for (auto j: endpoints_to_erase) {
            m_active_endpoints[key].erase(m_active_endpoints[key].begin() + j);
        }
    }

    void Constraintpoints::update_endpoints(Domain* domain) {
        remove_active_endpoint(domain);

        // Copy inactive endpoints to active endpoints
        bool endpoint_present {m_inactive_endpoints.count(domain) > 0};
        if (endpoint_present) {
            Domain* endpoint_domain {m_inactive_endpoints[domain]};
            add_active_endpoint(endpoint_domain, domain->m_pos);
        }
    }

    Domain* Constraintpoints::get_domain_to_grow(Domain* domain) {
        return m_growthpoints[domain];
    }

    bool Constraintpoints::endpoint_reached(Domain* domain, VectorThree pos) {
        bool reached {false};
        auto seg = m_segs[domain];
        pair<int, int> key {domain->m_c, seg};
        for (auto endpoint: m_active_endpoints[key]) {
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
            int dir,
            int offset) {

        long double prod_nws {1}; // Product of num ideal walks

        // Loop through all active endpoints on current chain
        pair<int, int> key {domain->m_c, m_segs[domain]};
        for (auto endpoint: m_active_endpoints[key]) {
            int end_d_i {endpoint.first};
            int steps {calc_remaining_steps(end_d_i, domain, dir, offset)};
            VectorThree end_p {endpoint.second};
            prod_nws *= m_ideal_random_walks.num_walks(pos, end_p, steps);
        }

        return prod_nws;
    }

    bool Constraintpoints::walks_remain(
            Domain* domain,
            VectorThree pos,
            int dir,
            int offset) {

        // Loop through all active endpoints on current chain
        bool w_remain {true};
        pair<int, int> key {domain->m_c, m_segs[domain]};
        for (auto endpoint: m_active_endpoints[key]) {
            int end_d_i {endpoint.first};
            int steps {calc_remaining_steps(end_d_i, domain, dir, offset)};
            VectorThree end_p {endpoint.second};
            if (m_ideal_random_walks.num_walks(pos, end_p, steps) == 0) {
                w_remain = false;
                break;
            }
        }

        return w_remain;
    }

    void Constraintpoints::find_growthpoints_endpoints(
            vector<Domain*> scaffold_domains,
            vector<int> excluded_staples,
            int seg) {
        
        for (auto d: scaffold_domains) {
            m_d_stack.push_back(d);
            m_segs[d] = seg;

            // Skip if not bound, bound to self or already checked
            bool bound {d->m_bound_domain != nullptr};
            if (not bound or
                    bound_to_self(d) or
                    chain_included(m_checked_staples, d->m_bound_domain->m_c)) {
                continue;
            }

            m_staple_network.scan_network(d, excluded_staples);
            auto net_cs = m_staple_network.get_participating_chains();
            auto pot_gps = m_staple_network.get_potential_growthpoints();
            auto pot_iaes = m_staple_network.get_potential_inactive_endpoints();
            auto pot_ds = m_staple_network.get_potential_domain_stack();
            auto s_seg_map = m_staple_network.get_staple_to_segs_map();
            if (m_staple_network.externally_bound()) {
                add_active_endpoints_on_scaffold(pot_gps, pot_iaes, seg);
            }
            else {
                add_growthpoints(pot_gps);
                add_inactive_endpoints(pot_iaes);
                add_regrowth_staples(net_cs, excluded_staples);
                add_domains_to_stack(pot_ds);
                add_staple_to_segs_maps(s_seg_map);
            }
            for (auto ci: net_cs) {
                m_checked_staples.insert(ci);
            }
        }
    }

    bool Constraintpoints::bound_to_self(Domain* d) {
        return d->m_bound_domain->m_c == d->m_c;
    }

    void Constraintpoints::add_growthpoints(
            vector<pair<Domain*, Domain*>> potential_growthpoints) {

        for (auto growthpoint: potential_growthpoints) {
            m_growthpoints[growthpoint.first] = growthpoint.second;
        }
    }

    void Constraintpoints::add_inactive_endpoints(
            vector<pair<Domain*, Domain*>> pot_iaes) {

        for (auto pot_iae: pot_iaes) {
            m_inactive_endpoints[pot_iae.first] = pot_iae.second;
        }
    }

    void Constraintpoints::add_regrowth_staples(
            set<int> participating_chains,
            vector<int> ex_staples) {

        for (auto c_i: participating_chains) {
            if (c_i == 0) {
                continue;
            }
            bool staple_excluded {chain_included(ex_staples, c_i)};
            if (staple_excluded) {
                m_regrowth_staples.insert(c_i);
            }
        }
    }

    void Constraintpoints::add_domains_to_stack(
            vector<Domain*> potential_d_stack) {

        m_d_stack.insert(m_d_stack.begin(), potential_d_stack.begin(),
                potential_d_stack.end());
    }

    void Constraintpoints::add_active_endpoints_on_scaffold(
            vector<pair<Domain*, Domain*>> pot_growthpoints,
            vector<pair<Domain*, Domain*>> pot_inactive_endpoints,
            int seg) {

        // Extract those from potential growthpoints
        for (auto growth_pair: pot_growthpoints) {
            Domain* growth_domain {growth_pair.first};
            if (growth_domain->m_c == 0) {
                add_active_endpoint(growth_domain, growth_domain->m_pos, seg);
            }
        }

        // Extract those from inactive endpoints
        for (auto pot_iae: pot_inactive_endpoints) {
            if (pot_iae.second->m_c == m_origami_system.c_scaffold) {
                add_active_endpoint(pot_iae.first, pot_iae.second->m_pos, seg);
            }
        }
    }

    void Constraintpoints::add_staple_to_segs_maps(
            unordered_map<Domain*, int> s_seg_map) {

        m_segs.insert(s_seg_map.begin(), s_seg_map.end());
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
}
