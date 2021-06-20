// cb_movetypes.cpp

#include <iostream>
#include <map>
#include <utility>

#include "algorithm"
#include "cb_movetypes.h"

namespace movetypes {

using std::map;

using utility::Occupancy;

CBMCMovetype::CBMCMovetype(
        OrigamiSystem& origami_system,
        RandomGens& random_gens,
        IdealRandomWalks& ideal_random_walks,
        vector<OrigamiOutputFile*> config_files,
        string label,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        MCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        RegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params) {}

void CBMCMovetype::reset_internal() {
    MCMovetype::reset_internal();
    m_bias = 1;
    m_new_bias = 1;
    m_new_modifier = 1;
    m_regrow_old = false;
}

void CBMCMovetype::add_external_bias() {
    m_ops.update_move_params();
    double total_bias {m_biases.calc_move()};
    m_bias *= exp(-total_bias);
}

void CBMCMovetype::calc_biases(
        const VectorThree p_prev,
        Domain& domain,
        configsT& configs,
        vector<double>& bfactors) {

    // Iterate through all possible new positions
    for (auto v: utility::vectors) {

        // Trial position vector
        VectorThree p_new {p_prev + v};

        // Check energies of each configuration
        switch (m_origami_system.position_occupancy(p_new)) {
        case Occupancy::bound:
        case Occupancy::misbound:
            continue;
        case Occupancy::unbound: {
            Domain* unbound_domain {m_origami_system.unbound_domain_at(p_new)};
            VectorThree o_new;
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

void CBMCMovetype::select_and_set_config(const int i, vector<Domain*> domains) {

    Domain* domain {domains[i]};
    Domain* prev_domain {domains[i - 1]};
    VectorThree p_prev {prev_domain->m_pos};

    // Calculate weights
    vector<pair<VectorThree, VectorThree>> configs {};
    vector<double> bfactors {};
    calc_biases(p_prev, *domain, configs, bfactors);
    vector<double> weights {calc_bias(bfactors, configs, domain)};
    if (m_rejected) {
        return;
    }

    if (not m_regrow_old) {
        select_and_set_new_config(weights, configs, *domain);
    }
    else {
        select_and_set_old_config(*domain);
    }
    write_config();

    // Reversion list update
    pair<int, int> key {domain->m_c, domain->m_d};
    m_assigned_domains.push_back(key);
}

void CBMCMovetype::select_and_set_new_config(
        const vector<double> weights,
        const configsT configs,
        Domain& domain) {

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

double CBMCMovetype::set_old_growth_point(
        Domain& growth_domain_new,
        Domain& growth_domain_old) {
    pair<int, int> key {growth_domain_new.m_c, growth_domain_new.m_d};
    VectorThree o_old {m_old_ore[key]};
    double delta_e {0};
    delta_e += m_origami_system.set_checked_domain_config(
            growth_domain_new, growth_domain_old.m_pos, o_old);
    m_bias *= exp(-delta_e);
    m_assigned_domains.push_back(key);

    return delta_e;
}

bool CBMCMovetype::test_cb_acceptance() {
    long double ratio {m_new_bias / m_bias};
    bool accepted;
    if (test_acceptance(ratio)) {
        reset_origami();
        m_ops.update_move_params();
        m_biases.calc_move();
        accepted = true;
    }
    else {
        m_modified_domains.clear();
        m_assigned_domains.clear();
        accepted = false;
    }

    return accepted;
}

double CBMCMovetype::unassign_domains(vector<Domain*> domains) {
    double delta_e {0};
    for (auto domain: domains) {
        pair<int, int> key {domain->m_c, domain->m_d};
        m_prev_pos[key] = domain->m_pos;
        m_prev_ore[key] = domain->m_ore;
        m_modified_domains.push_back(key);

        delta_e += m_origami_system.unassign_domain(*domain);
    }
    write_config();

    return delta_e;
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

CBStapleRegrowthMCMovetype::CBStapleRegrowthMCMovetype(
        OrigamiSystem& origami_system,
        RandomGens& random_gens,
        IdealRandomWalks& ideal_random_walks,
        vector<OrigamiOutputFile*> config_files,
        string label,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        MCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        RegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        CBMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params) {}

void CBStapleRegrowthMCMovetype::write_log_summary(ostream* log_stream) {
    write_log_summary_header(log_stream);

    // Insertion of each staple type
    map<int, int> attempts {};
    map<int, int> accepts {};
    set<int> staple_types {};
    for (auto tracker: m_tracking) {
        auto info = tracker.first;
        auto counts = tracker.second;
        if (not info.no_staples) {
            staple_types.insert(info.staple_type);
            attempts[info.staple_type] = counts.attempts;
            accepts[info.staple_type] = counts.accepts;
        }
    }
    *log_stream << "\n";
    for (auto st: staple_types) {
        *log_stream << "    Staple type: " << st << "\n";
        int ats {attempts[st]};
        int acs {accepts[st]};
        double freq {static_cast<double>(acs) / ats};
        *log_stream << "        Attempts: " << ats << "\n";
        *log_stream << "        Accepts: " << acs << "\n";
        *log_stream << "        Frequency: " << freq << "\n";
    }
}

bool CBStapleRegrowthMCMovetype::internal_attempt_move() {
    bool accepted {false};

    // No staples to regrow
    if (m_origami_system.num_staples() == 0) {
        accepted = false;
        m_tracker.no_staples = true;
        return accepted;
    }
    else {
        m_tracker.no_staples = false;
    }

    // Select a staple to regrow
    int c_i_index {
            m_random_gens.uniform_int(1, m_origami_system.num_staples())};
    vector<Domain*> selected_chain {m_origami_system.get_chains()[c_i_index]};
    m_tracker.staple_type = selected_chain[0]->m_c_ident;

    // Reject if staple is connector
    if (staple_is_connector(selected_chain)) {
        return accepted;
    }

    // Select growth points on chains
    auto bound_domains = find_bound_domains(selected_chain);
    m_bias *= bound_domains.size();
    domainPairT growthpoint {select_old_growthpoint(bound_domains)};
    unassign_domains(selected_chain);

    // Grow staple
    set_growthpoint_and_grow_staple(growthpoint, selected_chain);
    if (m_rejected) {
        return accepted;
    }
    m_bias /= num_bound_staple_domains(selected_chain);

    add_external_bias();

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

void CBStapleRegrowthMCMovetype::add_tracker(bool accepted) {
    movetypes::add_tracker(m_tracker, m_tracking, accepted);
}

void CBStapleRegrowthMCMovetype::grow_chain(vector<Domain*> domains) {
    for (size_t i {1}; i != domains.size(); i++) {
        select_and_set_config(i, domains);
        if (m_rejected) {
            break;
        }
    }
}

vector<double> CBStapleRegrowthMCMovetype::calc_bias(
        const vector<double> bfactors,
        const configsT&,
        Domain*) {

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

void CBStapleRegrowthMCMovetype::set_growthpoint_and_grow_staple(
        domainPairT growthpoint,
        vector<Domain*> selected_chain) {

    if (m_regrow_old) {
        set_old_growth_point(*growthpoint.first, *growthpoint.second);
    }
    else {
        double delta_e {
                set_growth_point(*growthpoint.first, *growthpoint.second)};
        m_bias *= exp(-delta_e);
    }
    if (not m_rejected) {
        grow_staple(growthpoint.first->m_d, selected_chain);
    }
}

CTCBRegrowthMCMovetype::CTCBRegrowthMCMovetype(
        OrigamiSystem& origami_system,
        RandomGens& random_gens,
        IdealRandomWalks& ideal_random_walks,
        vector<OrigamiOutputFile*> config_files,
        string label,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params,
        int num_excluded_staples,
        int max_regrowth):
        MCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        RegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        CBMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        CTRegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params,
                num_excluded_staples,
                max_regrowth,
                max_regrowth) {}

void CTCBRegrowthMCMovetype::reset_internal() {
    m_regrow_old = false;
    CBMCMovetype::reset_internal();
    CTRegrowthMCMovetype::reset_internal();
}

void CTCBRegrowthMCMovetype::grow_chain(vector<Domain*> domains) {
    if (domains.size() <= 1) {
        return;
    }
    for (size_t i {1}; i != domains.size(); i++) {
        Domain* domain {domains[i]};
        m_dir = m_constraintpoints.get_dir(domain);
        select_and_set_config(i, domains);
        if (m_rejected) {
            break;
        }
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

vector<double> CTCBRegrowthMCMovetype::calc_bias(
        const vector<double> bfactors,
        const vector<pair<VectorThree, VectorThree>>& configs,
        Domain* domain) {

    vector<long double> weights(bfactors.begin(), bfactors.end());
    for (size_t i {0}; i != configs.size(); i++) {
        VectorThree cur_pos {configs[i].first};

        // Set weights of non-self binding to 0 unless endpoint reached
        Occupancy pos_occ {m_origami_system.position_occupancy(cur_pos)};
        if (pos_occ == Occupancy::unbound) {
            Domain* occ_domain {m_origami_system.unbound_domain_at(cur_pos)};
            bool binding_same_chain {occ_domain->m_c == domain->m_c};
            bool endpoint {
                    m_constraintpoints.endpoint_reached(domain, cur_pos)};
            bool excluded_staple {find(m_excluded_staples.begin(),
                                       m_excluded_staples.end(),
                                       occ_domain->m_c) !=
                                  m_excluded_staples.end()};
            if (not(binding_same_chain or endpoint or excluded_staple)) {
                weights[i] = 0;
            }
        }

        // Zero probability of accepting configs with no walks remaining
        if (not m_constraintpoints.walks_remain(domain, cur_pos)) {
            weights[i] = 0;
        }
    }

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
        m_bias *= weights_sum;

        // Normalize
        for (size_t i {0}; i != weights.size(); i++) {
            double norm_weight = static_cast<double>(weights[i] / weights_sum);
            norm_weights.push_back(norm_weight);
        }
    }

    return norm_weights;
}

void CTCBRegrowthMCMovetype::grow_staple_and_update_endpoints(
        Domain* growth_d_old) {

    Domain* growth_d_new {m_constraintpoints.get_domain_to_grow(growth_d_old)};
    int c_i {growth_d_new->m_c};
    // HACK TO PREVENT ACCIDENTLY GROWING SCAFFOLD SEGMENTS
    if (c_i == m_origami_system.c_scaffold) {
        return;
    }
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

CTCBScaffoldRegrowthMCMovetype::CTCBScaffoldRegrowthMCMovetype(
        OrigamiSystem& origami_system,
        RandomGens& random_gens,
        IdealRandomWalks& ideal_random_walks,
        vector<OrigamiOutputFile*> config_files,
        string label,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params,
        int num_excluded_staples,
        int max_regrowth):
        MCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        RegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        CBMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        CTRegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params,
                num_excluded_staples,
                max_regrowth,
                max_regrowth),
        CTCBRegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params,
                num_excluded_staples,
                max_regrowth) {}

void CTCBScaffoldRegrowthMCMovetype::write_log_summary(ostream* log_stream) {
    write_log_summary_header(log_stream);
    // Insertion of each staple type
    map<int, int> length_attempts {};
    map<int, int> length_accepts {};
    map<int, int> staple_attempts {};
    map<int, int> staple_accepts {};
    set<int> lengths {};
    set<int> staples {};
    for (auto tracker: m_tracking) {
        auto info = tracker.first;
        auto counts = tracker.second;
        int length {info.num_scaffold_domains};
        int num_staples {info.num_staples};
        lengths.insert(length);
        staples.insert(num_staples);
        if (lengths.find(length) == lengths.end()) {
            length_attempts[length] = counts.attempts;
            length_accepts[length] = counts.accepts;
        }
        else {
            length_attempts[length] += counts.attempts;
            length_accepts[length] += counts.accepts;
        }
        if (staples.find(num_staples) == lengths.end()) {
            staple_attempts[num_staples] = counts.attempts;
            staple_accepts[num_staples] = counts.accepts;
        }
        else {
            staple_attempts[num_staples] += counts.attempts;
            staple_accepts[num_staples] += counts.accepts;
        }
    }
    for (auto l: lengths) {
        *log_stream << "    Number of scaffold domains: " << l << "\n";
        int ats {length_attempts[l]};
        int acs {length_accepts[l]};
        double freq {static_cast<double>(acs) / ats};
        *log_stream << "        Attempts: " << ats << "\n";
        *log_stream << "        Accepts: " << acs << "\n";
        *log_stream << "        Frequency: " << freq << "\n";
    }
    *log_stream << "\n";
    for (auto st: staples) {
        *log_stream << "    Number of staples: " << st << "\n";
        int ats {staple_attempts[st]};
        int acs {staple_accepts[st]};
        double freq {static_cast<double>(acs) / ats};
        *log_stream << "        Attempts: " << ats << "\n";
        *log_stream << "        Accepts: " << acs << "\n";
        *log_stream << "        Frequency: " << freq << "\n";
    }
}

bool CTCBScaffoldRegrowthMCMovetype::internal_attempt_move() {
    bool accepted {false};

    m_regrow_old = false;
    vector<Domain*> scaffold_domains {select_indices(m_scaffold, 2)};
    m_tracker.num_scaffold_domains = scaffold_domains.size();
    sel_excluded_staples();

    m_constraintpoints.calculate_constraintpoints(
            scaffold_domains, m_dir, m_excluded_staples);
    if (not(m_origami_system.m_cyclic and
            scaffold_domains.size() == m_origami_system.get_chain(0).size())) {
        m_constraintpoints.remove_active_endpoint(scaffold_domains[0]);
    }
    set<int> staples {m_constraintpoints.staples_to_be_regrown()};
    m_tracker.num_staples = staples.size();

    // Unassign domains to be regrown
    vector<Domain*> scaffold_domains_to_unassign(
            scaffold_domains.begin() + 1, scaffold_domains.end());
    unassign_domains(scaffold_domains_to_unassign);
    for (auto c_i: staples) {
        unassign_domains(m_origami_system.get_chain(c_i));
    }

    // Grow scaffold and staples
    if (m_constraintpoints.is_growthpoint(scaffold_domains[0])) {
        grow_staple_and_update_endpoints(scaffold_domains[0]);
        if (m_rejected) {
            return accepted;
        }
    }
    grow_chain(scaffold_domains);

    // Check if excluded staples have become unbound
    bool bound_to_system {excluded_staples_bound()};
    if (m_rejected or (m_num_excluded_staples != 0 and not bound_to_system)) {
        return accepted;
    }
    add_external_bias();

    // Regrow in old conformation
    setup_for_regrow_old();
    m_constraintpoints.reset_active_endpoints();
    if (not(m_origami_system.m_cyclic and
            scaffold_domains.size() == m_origami_system.get_chain(0).size())) {
        m_constraintpoints.remove_active_endpoint(scaffold_domains[0]);
    }

    // Unassign staples except those bound/linked to external scaffold
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

void CTCBScaffoldRegrowthMCMovetype::add_tracker(bool accepted) {
    movetypes::add_tracker(m_tracker, m_tracking, accepted);
}

CTCBJumpScaffoldRegrowthMCMovetype::CTCBJumpScaffoldRegrowthMCMovetype(
        OrigamiSystem& origami_system,
        RandomGens& random_gens,
        IdealRandomWalks& ideal_random_walks,
        vector<OrigamiOutputFile*> config_files,
        string label,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params,
        int num_excluded_staples,
        int max_regrowth,
        int max_seg_regrowth):
        MCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        RegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        CBMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params),
        CTRegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params,
                num_excluded_staples,
                max_regrowth,
                max_seg_regrowth),
        CTCBRegrowthMCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                config_files,
                label,
                ops,
                biases,
                params,
                num_excluded_staples,
                max_regrowth) {}

void CTCBJumpScaffoldRegrowthMCMovetype::write_log_summary(
        ostream* log_stream) {

    write_log_summary_header(log_stream);
    // Insertion of each staple type
    map<int, int> length_attempts {};
    map<int, int> length_accepts {};
    map<int, int> staple_attempts {};
    map<int, int> staple_accepts {};
    set<int> lengths {};
    set<int> staples {};
    for (auto tracker: m_tracking) {
        auto info = tracker.first;
        auto counts = tracker.second;
        int length {info.num_scaffold_domains};
        int num_staples {info.num_staples};
        lengths.insert(length);
        staples.insert(num_staples);
        if (lengths.find(length) == lengths.end()) {
            length_attempts[length] = counts.attempts;
            length_accepts[length] = counts.accepts;
        }
        else {
            length_attempts[length] += counts.attempts;
            length_accepts[length] += counts.accepts;
        }
        if (staples.find(num_staples) == lengths.end()) {
            staple_attempts[num_staples] = counts.attempts;
            staple_accepts[num_staples] = counts.accepts;
        }
        else {
            staple_attempts[num_staples] += counts.attempts;
            staple_accepts[num_staples] += counts.accepts;
        }
    }
    for (auto l: lengths) {
        *log_stream << "    Number of scaffold domains: " << l << "\n";
        int ats {length_attempts[l]};
        int acs {length_accepts[l]};
        double freq {static_cast<double>(acs) / ats};
        *log_stream << "        Attempts: " << ats << "\n";
        *log_stream << "        Accepts: " << acs << "\n";
        *log_stream << "        Frequency: " << freq << "\n";
    }
    *log_stream << "\n";
    for (auto st: staples) {
        *log_stream << "    Number of staples: " << st << "\n";
        int ats {staple_attempts[st]};
        int acs {staple_accepts[st]};
        double freq {static_cast<double>(acs) / ats};
        *log_stream << "        Attempts: " << ats << "\n";
        *log_stream << "        Accepts: " << acs << "\n";
        *log_stream << "        Frequency: " << freq << "\n";
    }
}

bool CTCBJumpScaffoldRegrowthMCMovetype::internal_attempt_move() {
    bool accepted {false};

    // Select scaffold segments and excluded staples
    m_regrow_old = false;
    vector<vector<Domain*>> segs {};
    vector<vector<vector<Domain*>>> paired_segs {};
    vector<Domain*> seg_stems {};
    vector<int> dirs {};
    select_noncontig_segs(m_scaffold, segs, paired_segs, seg_stems, dirs);
    m_tracker.num_scaffold_domains = 0;
    sel_excluded_staples();

    m_constraintpoints.calculate_constraintpoints(
            segs, dirs, m_excluded_staples);
    m_constraintpoints.remove_active_endpoint(segs[0][0]);
    set<int> staples {m_constraintpoints.staples_to_be_regrown()};
    m_tracker.num_staples = staples.size();
    vector<Domain*> regrow_domains {m_constraintpoints.domains_to_be_regrown()};

    // Unassign domains to be regrown
    vector<Domain*> domains_to_unassign(
            regrow_domains.begin() + 1, regrow_domains.end());
    unassign_domains(domains_to_unassign);

    // Grow scaffold and staples
    if (m_constraintpoints.is_growthpoint(segs[0][0])) {
        grow_staple_and_update_endpoints(segs[0][0]);
        if (m_rejected) {
            return accepted;
        }
    }
    grow_chain(segs[0]);
    if (m_rejected) {
        return accepted;
    }
    for (size_t i {0}; i != seg_stems.size(); i++) {
        Domain* seg_stem {seg_stems[i]};
        set_first_seg_domain(seg_stem);
        if (m_rejected) {
            return accepted;
        }
        for (auto seg: paired_segs[i]) {
            seg.insert(seg.begin(), seg_stem);
            grow_chain(seg);
            if (m_rejected) {
                return accepted;
            }
        }
    }

    // Check if excluded staples have become unbound
    bool bound_to_system {excluded_staples_bound()};
    if (m_rejected or (m_num_excluded_staples != 0 and not bound_to_system)) {
        return accepted;
    }
    add_external_bias();

    // Regrow in old conformation
    setup_for_regrow_old();
    m_constraintpoints.reset_active_endpoints();
    m_constraintpoints.remove_active_endpoint(segs[0][0]);
    unassign_domains(domains_to_unassign);

    // Grow scaffold and staples
    if (m_constraintpoints.is_growthpoint(segs[0][0])) {
        grow_staple_and_update_endpoints(segs[0][0]);
        if (m_rejected) {
            return accepted;
        }
    }
    grow_chain(segs[0]);
    if (m_rejected) {
        return accepted;
    }
    for (size_t i {0}; i != seg_stems.size(); i++) {
        Domain* seg_stem {seg_stems[i]};
        set_first_seg_domain(seg_stem);
        if (m_rejected) {
            return accepted;
        }
        for (auto seg: paired_segs[i]) {
            seg.insert(seg.begin(), seg_stem);
            grow_chain(seg);
            if (m_rejected) {
                return accepted;
            }
        }
    }

    // Reset modifier and test acceptance
    m_modifier = 1;
    accepted = test_cb_acceptance();

    return accepted;
}

void CTCBJumpScaffoldRegrowthMCMovetype::add_tracker(bool accepted) {
    movetypes::add_tracker(m_tracker, m_tracking, accepted);
}

void CTCBJumpScaffoldRegrowthMCMovetype::set_first_seg_domain(Domain* d_new) {

    Domain* d_old {m_constraintpoints.get_growthpoint(d_new)};
    if (not m_constraintpoints.walks_remain(d_new, d_old->m_pos)) {
        m_rejected = true;
        return;
    }
    m_dir = m_constraintpoints.get_dir(d_new);
    double delta_e {set_growth_point(*d_new, *d_old)};
    m_bias *= exp(-delta_e);
}
} // namespace movetypes
