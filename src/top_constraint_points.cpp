// top_constraint_points.cpp

#include <cmath>

#include "top_constraint_points.h"
#include "utility.h"

namespace topConstraintPoints {

using std::cout;
using std::find;
using std::set;

bool domain_included(vector<Domain*> domains, Domain* d) {
    return find(domains.begin(), domains.end(), d) != domains.end();
}

bool chain_included(vector<int> staples, int staple) {
    return find(staples.begin(), staples.end(), staple) != staples.end();
}

bool chain_included(set<int> staples, int staple) {
    return staples.find(staple) != staples.end();
}

StapleNetwork::StapleNetwork(OrigamiSystem& origami): m_origami {origami} {}

void StapleNetwork::set_excluded_staples(vector<int> excluded_staples) {
    m_ex_staples = excluded_staples;
}

void StapleNetwork::set_scaffold_domains(vector<Domain*> scaffold_domains) {
    m_scaffold_ds = scaffold_domains;
}

void StapleNetwork::scan_network(Domain* d) {
    clear_network();

    // Adding the scaffold chain made checks easier
    m_net_cs.insert(m_origami.c_scaffold);

    // Start recursive scan of network
    scan_staple_topology(d);
}

set<int> StapleNetwork::get_participating_chains() { return m_net_cs; }

vector<pair<Domain*, Domain*>> StapleNetwork::get_potential_growthpoints() {
    return m_pot_gps;
}

vector<pair<Domain*, Domain*>> StapleNetwork::
        get_potential_inactive_endpoints() {
    return m_pot_iaes;
}

vector<Domain*> StapleNetwork::get_potential_domain_stack() { return m_pot_ds; }

unordered_map<Domain*, int> StapleNetwork::get_staple_to_segs_map() {
    return m_segs;
}

unordered_map<pair<int, int>, int> StapleNetwork::get_dirs() {
    return m_domain_to_dir;
}

bool StapleNetwork::externally_bound() { return m_external; }

void StapleNetwork::clear_network() {
    m_external = false;
    m_net_cs.clear();
    m_pot_gps.clear();
    m_pot_ds.clear();
    m_pot_iaes.clear();
    m_segs.clear();
}

void StapleNetwork::scan_staple_topology(Domain* growth_d) {

    int ci {growth_d->m_c}; // Staple chain index
    m_net_cs.insert(ci);
    m_pot_ds.push_back(growth_d);
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

        /* If the staple has not been added to the network, recurse into
         * the bound staple. Otherwise consider adding potential endpoints
         */
        bool bs_excluded {chain_included(m_ex_staples, bd_ci)};
        bool in_network {chain_included(m_net_cs, bd_ci)};
        if (in_network) {

            /* If the bound domain is on the scaffold and in the range
             * being regrown or on another staple, and nothing is excluded
             * the domans are potential endpoints. Otherwise the network is
             * externally bound
             */
            bool d_excluded {false};
            bool external {false};
            if (bd_ci == 0) {
                external = not domain_included(m_scaffold_ds, bd);
            }
            else {
                d_excluded = chain_included(m_ex_staples, ci);
            }
            if (not external and not(bs_excluded and d_excluded)) {
                add_potential_inactive_endpoint(d, bd);
            }
            if (not m_external and external) {
                m_external = external;
            }
        }
        else {

            /* Add domain to potential growthpoints if bound domain if not
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
    pair<int, int> key {ci, seg};
    m_domain_to_dir[key] = 1;
    for (size_t di = d->m_d + 1; di != staple.size(); di++) {
        Domain* d_loop {staple[di]};
        ds.push_back(d_loop);
        m_segs[d_loop] = seg;
    }

    // Add domains in five prime direction
    seg++;
    key = {ci, seg};
    m_domain_to_dir[key] = -1;
    for (int di {d->m_d - 1}; di != -1; di--) {
        Domain* d_loop {staple[di]};
        ds.push_back(d_loop);
        m_segs[d_loop] = seg;
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
        IdealRandomWalks& ideal_random_walks):
        m_origami_system {origami_system},
        m_ideal_random_walks {ideal_random_walks},
        m_staple_network {origami_system} {}

int Constraintpoints::get_dir(Domain* d) {
    auto seg = m_segs[d];
    pair<int, int> key {d->m_c, seg};
    int dir {m_domain_to_dir[key]};

    return dir;
}

vector<VectorThree> Constraintpoints::get_erased_endpoints() {
    return m_erased_endpoints;
}

void Constraintpoints::reset_internal() {
    m_scaffold_domains.clear();
    m_regrowth_staples.clear();
    m_checked_staples.clear();
    m_erased_endpoints.clear();
    m_d_stack.clear();
    m_segs.clear();
    m_domain_to_dir.clear();
    m_growthpoints.clear();
    m_stemdomains.clear();
    m_active_endpoints.clear();
    m_initial_active_endpoints.clear();
    m_inactive_endpoints.clear();
    m_stemd_to_segs.clear();
}

void Constraintpoints::calculate_constraintpoints(
        vector<Domain*> scaffold_domains,
        int dir,
        vector<int> excluded_staples) {

    m_scaffold_domains = scaffold_domains;
    m_staple_network.set_scaffold_domains(m_scaffold_domains);
    m_staple_network.set_excluded_staples(excluded_staples);
    pair<int, int> key {scaffold_domains[0]->m_c, 0};
    m_domain_to_dir[key] = dir;
    find_growthpoints_endpoints(scaffold_domains, excluded_staples, 0);

    // Save initial active endpoints and remaining steps for regrowing old
    m_initial_active_endpoints = m_active_endpoints;
}

void Constraintpoints::calculate_constraintpoints(
        vector<vector<Domain*>> scaffold_segments,
        vector<int> dirs,
        vector<int> excluded_staples) {

    for (auto segment: scaffold_segments) {
        m_scaffold_domains.insert(
                m_scaffold_domains.end(), segment.begin(), segment.end());
    }
    m_staple_network.set_scaffold_domains(m_scaffold_domains);
    m_staple_network.set_excluded_staples(excluded_staples);
    int seg {0};
    for (size_t i {0}; i != scaffold_segments.size(); i++) {
        auto scaffold_domains = scaffold_segments[i];
        if (scaffold_domains.size() != 0) {
            pair<int, int> key {scaffold_domains[0]->m_c, seg};
            m_domain_to_dir[key] = dirs[i];
            find_growthpoints_endpoints(
                    scaffold_domains, excluded_staples, seg);
        }
        seg++;
    }

    // Save initial active endpionts and remaining steps for regrowing old
    m_initial_active_endpoints = m_active_endpoints;
}

set<int> Constraintpoints::staples_to_be_regrown() {
    return m_regrowth_staples;
}

vector<Domain*> Constraintpoints::domains_to_be_regrown() { return m_d_stack; }

bool Constraintpoints::is_growthpoint(Domain* domain) {
    return (m_growthpoints.count(domain) > 0);
}

bool Constraintpoints::is_stemdomain(Domain* domain) {
    return (m_stemdomains.count(domain) > 0);
}

void Constraintpoints::add_active_endpoint(Domain* d, VectorThree pos) {

    auto seg = m_segs.at(d);
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

void Constraintpoints::add_inactive_endpoint(Domain* d_i, Domain* d_j) {
    m_inactive_endpoints[d_i] = d_j;
}

void Constraintpoints::add_growthpoint(Domain* growthpoint, Domain* stemd) {

    m_growthpoints[growthpoint] = stemd;
    m_stemdomains[stemd] = growthpoint;
}

void Constraintpoints::add_stem_seg_pair(Domain* stemd, vector<int> seg_pair) {

    m_stemd_to_segs[stemd] = seg_pair;
}

void Constraintpoints::reset_active_endpoints() {
    m_active_endpoints = m_initial_active_endpoints;
}

void Constraintpoints::remove_active_endpoint(Domain* domain) {

    // Remove endpoint if reached
    m_erased_endpoints.clear();
    int c_i {domain->m_c};
    auto seg = m_segs.at(domain);
    pair<int, int> key {c_i, seg};
    vector<int> endpoints_to_erase {};
    for (size_t j {0}; j != m_active_endpoints[key].size(); j++) {
        pair<int, VectorThree> endpoint {m_active_endpoints[key][j]};
        if (endpoint.first == domain->m_d) {
            endpoints_to_erase.push_back(j);
            m_erased_endpoints.push_back(endpoint.second);
        }
    }

    // Removing in reverse means the indices are not changing as I do it
    std::reverse(endpoints_to_erase.begin(), endpoints_to_erase.end());
    for (auto j: endpoints_to_erase) {
        m_active_endpoints[key].erase(m_active_endpoints[key].begin() + j);
    }
}

void Constraintpoints::remove_activated_endpoint(Domain* domain) {
    bool endpoint_present {m_inactive_endpoints.count(domain) > 0};
    if (endpoint_present) {
        Domain* edomain {m_inactive_endpoints[domain]};
        int c_i {edomain->m_c};
        auto seg = m_segs.at(edomain);
        pair<int, int> key {c_i, seg};
        for (size_t j {0}; j != m_active_endpoints[key].size(); j++) {
            pair<int, VectorThree> endpoint {m_active_endpoints[key][j]};
            if (endpoint.first == edomain->m_d) {
                m_active_endpoints[key].erase(
                        m_active_endpoints[key].begin() + j);
                break;
            }
        }
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

Domain* Constraintpoints::get_growthpoint(Domain* domain) {
    return m_stemdomains[domain];
}

bool Constraintpoints::endpoint_reached(Domain* domain, VectorThree pos) {
    bool reached {false};
    auto seg = m_segs.at(domain);
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

    // Loop through all active endpoints on current segment
    pair<int, int> key {domain->m_c, m_segs.at(domain)};
    for (auto endpoint: m_active_endpoints[key]) {
        int end_d_i {endpoint.first};
        int steps {calc_remaining_steps(end_d_i, domain, dir, offset)};
        if (steps < 0) {
            cout << "Bad endpoint detected\n";
            throw utility::SimulationMisuse {};
        }
        VectorThree end_p {endpoint.second};
        prod_nws *= m_ideal_random_walks.num_walks(pos, end_p, steps);
    }

    return prod_nws;
}

bool Constraintpoints::walks_remain(
        Domain* domain,
        VectorThree pos,
        int offset) {

    vector<int> segs {};
    bool w_remain {true};
    if (is_stemdomain(domain)) {
        segs = m_stemd_to_segs[domain];
    }
    else {
        segs = {m_segs.at(domain)};
    }
    for (auto seg: segs) {
        pair<int, int> key {domain->m_c, seg};
        int dir {m_domain_to_dir[key]};
        w_remain = walks_remain(key, domain, pos, dir, offset);
        if (not w_remain) {
            break;
        }
    }

    return w_remain;
}

bool Constraintpoints::walks_remain(
        pair<int, int> key,
        Domain* domain,
        VectorThree pos,
        int dir,
        int offset) {

    // Loop through all active endpoints on current segment
    bool w_remain {true};
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

vector<pair<int, VectorThree>> Constraintpoints::get_active_endpoints(
        int c_i,
        int seg) {
    return m_active_endpoints[{c_i, seg}];
}

Domain* Constraintpoints::get_inactive_endpoints(Domain* domain) {
    return m_inactive_endpoints[domain];
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
        if (not bound or bound_to_self(d) or
            chain_included(m_checked_staples, d->m_bound_domain->m_c)) {
            continue;
        }

        Domain* bd {d->m_bound_domain};
        m_staple_network.scan_network(bd);
        auto net_cs = m_staple_network.get_participating_chains();
        auto pot_gps = m_staple_network.get_potential_growthpoints();
        pot_gps.push_back({d, bd});
        auto pot_iaes = m_staple_network.get_potential_inactive_endpoints();
        auto pot_ds = m_staple_network.get_potential_domain_stack();
        auto s_seg_map = m_staple_network.get_staple_to_segs_map();
        auto domain_to_dirs = m_staple_network.get_dirs();
        if (m_staple_network.externally_bound()) {
            add_active_endpoints_on_scaffold(pot_gps, pot_iaes, seg);
        }
        else {
            add_growthpoints(pot_gps);
            add_inactive_endpoints(pot_iaes);
            add_regrowth_staples(net_cs, excluded_staples);
            add_domains_to_stack(pot_ds);
            add_staple_to_segs_maps(s_seg_map);
            add_dirs(domain_to_dirs);
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
        m_stemdomains[growthpoint.second] = growthpoint.first;
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
        if (not staple_excluded) {
            m_regrowth_staples.insert(c_i);
        }
    }
}

void Constraintpoints::add_domains_to_stack(vector<Domain*> potential_d_stack) {

    m_d_stack.insert(
            m_d_stack.end(),
            potential_d_stack.begin(),
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

void Constraintpoints::add_dirs(unordered_map<pair<int, int>, int> dirs) {
    m_domain_to_dir.insert(dirs.begin(), dirs.end());
}

int Constraintpoints::calc_remaining_steps(
        int endpoint_d_i,
        Domain* domain,
        int dir,
        int step_offset) {
    int steps;
    if (m_origami_system.m_cyclic and
        domain->m_c == m_origami_system.c_scaffold) {

        if (dir > 0 and endpoint_d_i < domain->m_d) {
            steps = m_origami_system.get_chain(0).size() + endpoint_d_i -
                    domain->m_d;
        }
        else if (dir < 0 and endpoint_d_i > domain->m_d) {
            steps = domain->m_d + m_origami_system.get_chain(0).size() -
                    endpoint_d_i;
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
} // namespace topConstraintPoints
