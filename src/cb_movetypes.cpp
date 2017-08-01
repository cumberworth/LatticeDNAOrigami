// cb_movetypes.cpp

#include <iostream>
#include <map>
#include <utility>

#include "algorithm"
#include "cb_movetypes.h"

namespace movetypes {

    using std::cout;
    using std::map;

    using utility::Occupancy;
    using utility::OrigamiMisuse;

    void CBMCMovetype::reset_internal() {
        MCMovetype::reset_internal();
        m_bias = 1;
        m_new_bias = 1;
        m_new_modifier = 1;
        m_regrow_old = false;
        m_old_pos.clear();
        m_old_ore.clear();
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

    void CBMCMovetype::select_and_set_config(
            const int i,
            vector<Domain*> domains) {

        Domain* domain {domains[i]};
        Domain* prev_domain {domains[i - 1]};
        VectorThree p_prev {prev_domain->m_pos};

        // Calculate weights
        vector<pair<VectorThree, VectorThree>> configs {};
        vector<double> bfactors {};
        calc_biases(p_prev, *domain, configs, bfactors);
        vector<double> weights {calc_bias(bfactors, configs, p_prev, domain, domains)};
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
        delta_e += m_origami_system.set_checked_domain_config(growth_domain_new,
                growth_domain_old.m_pos, o_old);
        m_bias *= exp(-delta_e);
        m_assigned_domains.push_back(key);

        return delta_e;
    }

    bool CBMCMovetype::test_cb_acceptance() {
        long double ratio {m_new_bias / m_bias};
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
        write_config();
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

    vector<domainPairT> CBMCMovetype::find_bound_domains(
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
            cout << "System has unbound staple\n";
            throw OrigamiMisuse {};
        }

        return bound_domains;
    }

    domainPairT CBMCMovetype::select_old_growthpoint(
            vector<domainPairT> bound_domains) {

        int bound_domain_index {m_random_gens.uniform_int(0, bound_domains.size() - 1)};
        Domain* growth_domain_new {bound_domains[bound_domain_index].first};
        Domain* growth_domain_old {bound_domains[bound_domain_index].second};
        return {growth_domain_new, growth_domain_old};
    }

    void CBMCMovetype::add_external_bias() {
        m_ops.update_move_params();
        double total_bias {m_biases.calc_move()};
        m_bias *= exp(-total_bias);
    }

    void CBMCMovetype::update_external_bias() {
        m_ops.update_move_params();
        m_biases.calc_move();
    }

    bool CBStapleRegrowthMCMovetype::attempt_move(long long int step) {
        m_step = step;
        write_config();
        bool accepted {false};
        m_general_tracker.attempts++;

        // No staples to regrow
        if (m_origami_system.num_staples() == 0) {
            accepted = false;
            m_tracker.no_staples = true;
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }
        else {
            m_tracker.no_staples = false;
        }

        // Select a staple to regrow
        int c_i_index {m_random_gens.uniform_int(1, m_origami_system.num_staples())};
        vector<Domain*> selected_chain {m_origami_system.get_chains()[c_i_index]};
        m_tracker.staple_type = selected_chain[0]->m_c_ident;

        // Reject if staple is connector
        if (staple_is_connector(selected_chain)) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        // Select growth points on chains
        pair<Domain*, Domain*> growthpoint {select_new_growthpoint(selected_chain)};

        auto bound_domains = find_bound_domains(selected_chain);
        unassign_domains(selected_chain);
        update_external_bias();

        // Grow staple
        set_growthpoint_and_grow_staple(growthpoint, selected_chain);
        if (m_rejected) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        add_external_bias();

        // Regrow staple in old conformation
        setup_for_regrow_old();

        // Select growth point from previously bound domains
        growthpoint = select_old_growthpoint(bound_domains);

        // Unassign and add to reversion list
        unassign_domains(selected_chain);
        update_external_bias();

        // Grow staple
        set_growthpoint_and_grow_staple(growthpoint, selected_chain);

        // Revert modifier and test acceptance
        m_modifier = m_new_modifier;
        add_external_bias();
        accepted = test_cb_acceptance();
        add_tracker(m_tracker, m_tracking, accepted);
        m_general_tracker.accepts += accepted;

        return accepted;
    }

    void CBStapleRegrowthMCMovetype::write_log_summary(ostream* log_stream) {

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

    void CBStapleRegrowthMCMovetype::set_growthpoint_and_grow_staple(
            domainPairT growthpoint,
            vector<Domain*> selected_chain) {

        if (m_regrow_old) {
            set_old_growth_point(*growthpoint.first, *growthpoint.second);
        }
        else {
            double delta_e {set_growth_point(*growthpoint.first, *growthpoint.second)};
            m_bias *= exp(-delta_e);
        }
        if (not m_rejected) {
            grow_staple(growthpoint.first->m_d, selected_chain);
        }
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
            const VectorThree,
            Domain*,
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

    CTCBRegrowthMCMovetype::CTCBRegrowthMCMovetype(
            OrigamiSystem& origami_system,
            RandomGens& random_gens,
            IdealRandomWalks& ideal_random_walks,
            vector<OrigamiOutputFile*> config_files,
            string label,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params,
            int num_excluded_staples):
            CBMCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            m_num_excluded_staples {num_excluded_staples} {
    }

    void CTCBRegrowthMCMovetype::reset_internal() {
        CBMCMovetype::reset_internal();
        m_constraintpoints.reset_internal();
        m_excluded_staples.clear();
    }

    pair<int, int> CTCBRegrowthMCMovetype::select_endpoints(
            const int array_size,
            const int m,
            const int min_size) {

        int start_i {m_random_gens.uniform_int(m, array_size - 1 - m)};
        int end_i {m_random_gens.uniform_int(m, array_size - 1 - m)};
        while (std::abs(start_i - end_i) < min_size) {
            end_i = m_random_gens.uniform_int(m, array_size - 1 - m);
        }

        return {start_i, end_i};
    }

    vector<Domain*> CTCBRegrowthMCMovetype::select_linear_segment(
            Domain* start_domain,
            Domain* end_domain) {

        vector<Domain*> scaffold {m_origami_system.get_chain(m_origami_system.c_scaffold)};
        vector<Domain*> domains {};

        // Find direction of regrowth
        if (end_domain->m_d > start_domain->m_d) {
            m_dir = 1;
        }
        else {
            m_dir = -1;
        }
        for (int d_i {start_domain->m_d}; d_i != end_domain->m_d + m_dir; d_i += m_dir) {
            Domain* cur_domain {scaffold[d_i]};
            domains.push_back(cur_domain);
        }

        return domains;
    }

    vector<Domain*> CTCBRegrowthMCMovetype::select_cyclic_segment(
            Domain* start_domain,
            Domain* end_domain) {

        vector<Domain*> domains {};

        // Select direction of regrowth MESSY
        m_dir = m_random_gens.uniform_int(0, 1);
        if (m_dir == 0) {
            m_dir = -1;
        }
        int d_i {start_domain->m_d};
        Domain* cur_domain {start_domain};
        while (d_i != end_domain->m_d) {
            domains.push_back(cur_domain);
            cur_domain = (*cur_domain) + m_dir;
            d_i = cur_domain->m_d;
        }
        domains.push_back(end_domain);

        return domains;
    }

    void CTCBRegrowthMCMovetype::grow_chain(vector<Domain*> domains) {
        m_dir = domains[1]->m_d - domains[0]->m_d;
        for (size_t i {1}; i != domains.size(); i++) {
            select_and_set_config(i, domains);
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

    void CTCBRegrowthMCMovetype::grow_staple_and_update_endpoints(
            Domain* growth_d_old) {

        Domain* growth_d_new {m_constraintpoints.get_domain_to_grow(growth_d_old)};
        int c_i {growth_d_new->m_c};
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

    vector<double> CTCBRegrowthMCMovetype::calc_bias(
            const vector<double> bfactors,
            const vector<pair<VectorThree, VectorThree>>& configs,
            const VectorThree p_prev,
            Domain* domain,
            vector<Domain*> domains) {


        vector<long double> weights(bfactors.begin(), bfactors.end());
        for (size_t i {0}; i != configs.size(); i++) {
            VectorThree cur_pos {configs[i].first};

            // Set weights of positions involving non-self binding to 0 unless endpoint reached
            Occupancy pos_occ {m_origami_system.position_occupancy(cur_pos)};
            if (pos_occ == Occupancy::unbound) {
                Domain* occ_domain {m_origami_system.unbound_domain_at(cur_pos)};
                bool binding_same_chain {occ_domain->m_c == domain->m_c};
                bool endpoint {m_constraintpoints.endpoint_reached(domain, cur_pos)};
                bool excluded_staple {find(m_excluded_staples.begin(),
                        m_excluded_staples.end(), occ_domain->m_c) !=
                        m_excluded_staples.end()};
                if (not (binding_same_chain or endpoint or excluded_staple)) {
                    weights[i] = 0;
                }
            }

            // Bias weights with number of walks
            weights[i] *= m_constraintpoints.calc_num_walks_prod(domain, cur_pos,
                    domains, m_dir);
        }

        // Calculate number of walks for previous position
        long double prod_num_walks {m_constraintpoints.calc_num_walks_prod(domain,
                p_prev, domains, m_dir, 1)};

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
            long double bias {weights_sum / prod_num_walks};
            m_bias *= bias;

            // Normalize
            for (size_t i {0}; i != weights.size(); i++) {
                double norm_weight = static_cast<double>(weights[i] / weights_sum);
                norm_weights.push_back(norm_weight);
            }
        }

        return norm_weights;
    }

    bool CTCBScaffoldRegrowthMCMovetype::attempt_move(long long int step) {
        m_step = step;
        write_config();
        bool accepted {false};
        m_general_tracker.attempts++;

        m_regrow_old = false;
        vector<Domain*> scaffold {m_origami_system.get_chain(
                m_origami_system.c_scaffold)};
        vector<Domain*> scaffold_domains {select_indices(scaffold)};
        m_tracker.num_scaffold_domains = scaffold_domains.size();

        for (int i {0}; i != m_num_excluded_staples; i++) {
            if (i == m_origami_system.num_staples()) {
                break;
            }
            bool staple_picked {false};
            int s_ui;
            while (not staple_picked)  {
                int staple_i {m_random_gens.uniform_int(1,
                        m_origami_system.num_staples())};
                s_ui = m_origami_system.get_chains()[staple_i][0]->m_c;
                staple_picked = (find(m_excluded_staples.begin(),
                        m_excluded_staples.end(), s_ui) ==
                        m_excluded_staples.end());
            }
            m_excluded_staples.push_back(s_ui);
        }
        m_constraintpoints.calculate_constraintpoints(scaffold_domains, m_excluded_staples);
        if (not m_origami_system.m_cyclic and scaffold_domains.size() !=
            m_origami_system.get_chain(0).size()) {
            m_constraintpoints.remove_active_endpoint(scaffold_domains[0]);
        }
        set<int> staples {m_constraintpoints.staples_to_be_regrown()};
        m_tracker.num_staples = staples.size();

        // Unassign staples except those directly and indirectly bound to external scaffold domains
        vector<Domain*> scaffold_domains_to_unassign(scaffold_domains.begin() + 1,
                scaffold_domains.end());
        unassign_domains(scaffold_domains_to_unassign);
        for (auto c_i: staples) {
            unassign_domains(m_origami_system.get_chain(c_i));
        }
        update_external_bias();

        // Grow scaffold and staples
        if (m_constraintpoints.is_growthpoint(scaffold_domains[0])) {
            grow_staple_and_update_endpoints(scaffold_domains[0]);
            if (m_rejected) {
                add_tracker(m_tracker, m_tracking, accepted);
                return accepted;
            }
        }
        grow_chain(scaffold_domains);

        // Check if excluded staples have become unbound
        bool bound_to_system {false};
        for (auto exs_i: m_excluded_staples) {
            vector<Domain*> exs {m_origami_system.get_chain(exs_i)};
            for (auto exd: exs) {
                set<int> dummy_set {};
                if (scan_for_scaffold_domain(exd, dummy_set)) {
                    bound_to_system = true;
                    break;
                }
            }
        }
        if (m_rejected or not bound_to_system) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        add_external_bias();
        // Regrow in old conformation
        setup_for_regrow_old();
        m_constraintpoints.reset_active_endpoints();
        if (not m_origami_system.m_cyclic and scaffold_domains.size() !=
            m_origami_system.get_chain(0).size()) {
            m_constraintpoints.remove_active_endpoint(scaffold_domains[0]);
        }

        // Unassign staples except those directly and indirectly bound to external scaffold domains
        unassign_domains(scaffold_domains_to_unassign);
        for (auto c_i: staples) {
            unassign_domains(m_origami_system.get_chain(c_i));
        }
        update_external_bias();

        // Grow scaffold and staples
        if (m_constraintpoints.is_growthpoint(scaffold_domains[0])) {
            grow_staple_and_update_endpoints(scaffold_domains[0]);
        }

        grow_chain(scaffold_domains);

        // Reset modifier and test acceptance
        m_modifier = 1;
        add_external_bias();
        accepted = test_cb_acceptance();
        add_tracker(m_tracker, m_tracking, accepted);
        m_general_tracker.accepts += accepted;
        return accepted;
    }

    void CTCBScaffoldRegrowthMCMovetype::write_log_summary(ostream* log_stream) {
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

    vector<Domain*> CTCBScaffoldRegrowthMCMovetype::select_indices(
            vector<Domain*> segment) {

        pair<int, int> endpoints {select_endpoints(segment.size(), 0, 1)};
        Domain* start_domain {segment[endpoints.first]};
        Domain* end_domain {segment[endpoints.second]};
        vector<Domain*> domains {};
        if (m_origami_system.m_cyclic) {
            domains = select_cyclic_segment(start_domain, end_domain);
        }
        else {
            domains = select_linear_segment(start_domain, end_domain);
        }

        // If end domain is end of chain, no endpoint
        Domain* endpoint_domain {(*end_domain) + m_dir};
        if (endpoint_domain != nullptr) {
            m_constraintpoints.add_active_endpoint(endpoint_domain, endpoint_domain->m_pos);
        }

        return domains;
    }

    CTCBLinkerRegrowthMCMovetype::CTCBLinkerRegrowthMCMovetype(
            OrigamiSystem& origami_system,
            RandomGens& random_gens,
            IdealRandomWalks& ideal_random_walks,
            vector<OrigamiOutputFile*> config_files,
            string label,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params,
            int num_excluded_staples,
            int max_disp,
            int max_turns):
            CTCBRegrowthMCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params, num_excluded_staples),
            m_max_disp {max_disp},
            m_max_turns {max_turns} {
    }

    bool CTCBLinkerRegrowthMCMovetype::attempt_move(long long int step) {
        m_step = step;
        write_config();
        bool accepted {false};
        m_regrow_old = false;
        m_general_tracker.attempts++;

        // Select region to be transformed and its linkers
        vector<Domain*> linker1 {};
        vector<Domain*> linker2 {};
        vector<Domain*> central_segment {};
        set<int> staples {select_and_setup_segments(linker1, linker2,
                central_segment)};
        set<int> central_staples {find_staples(central_segment)};
        vector<Domain*> central_domains {central_segment};
        for (auto staple: central_staples) {
            vector<Domain*> staple_domains {m_origami_system.get_chain(staple)};
            for (auto domain: staple_domains) {
                central_domains.push_back(domain);
            }
        }
        m_tracker.num_linker_domains = linker1.size() + linker2.size();
        m_tracker.num_linker_staples = staples.size();
        m_tracker.num_central_domains = central_segment.size();
        m_tracker.num_central_staples = central_staples.size();

        // Reject moves that have central region bound externally
        if (domains_bound_externally(central_segment)) {
            m_tracker.central_domains_connected = true;
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        // Unassign domains
        vector<Domain*> linker1_to_unassign(linker1.begin() + 1, linker1.end());
        unassign_domains(linker1_to_unassign);
        vector<Domain*> linker2_to_unassign(linker2.begin() + 1, linker2.end());
        unassign_domains(linker2_to_unassign);
        for (auto c_i: staples) {
            unassign_domains(m_origami_system.get_chain(c_i));
        }
        unassign_domains(central_domains); // Consider a more efficient method for this
        update_external_bias();

        transform_segment(linker1, linker2, central_segment, central_domains);

        // Grow linkers
        for (auto linker: {linker1, linker2}) {
            if (m_constraintpoints.is_growthpoint(linker[0])) {
                grow_staple_and_update_endpoints(linker[0]);
                if (m_rejected) {
                    add_tracker(m_tracker, m_tracking, accepted);
                    return accepted;
                }
            }
            grow_chain(linker);
            if (m_rejected) {
                add_tracker(m_tracker, m_tracking, accepted);
                return accepted;
            }
        }

        // Regrow in old conformation
        add_external_bias();
        setup_for_regrow_old();
        m_constraintpoints.reset_active_endpoints();

        // Unassign domains
        unassign_domains(linker1_to_unassign);
        unassign_domains(linker2_to_unassign);
        for (auto c_i: staples) {
            unassign_domains(m_origami_system.get_chain(c_i));
        }
        unassign_domains(central_domains); // Consider a more efficient method for this
        update_external_bias();

        revert_transformation(central_domains);

        // Grow linkers
        for (auto linker: {linker1, linker2}) {
            if (m_constraintpoints.is_growthpoint(linker[0])) {
                grow_staple_and_update_endpoints(linker[0]);
            }
            grow_chain(linker);
        }

        // Reset modifier and test acceptance
        m_modifier = 1;
        add_external_bias();
        accepted = test_cb_acceptance();
        add_tracker(m_tracker, m_tracking, accepted);
        m_general_tracker.accepts += accepted;
        return accepted;
    }

    void CTCBLinkerRegrowthMCMovetype::write_log_summary(ostream* log_stream) {
        // Insertion of each staple type
        map<pair<int, int>, pair<int, int>> length_counts {};
        set<pair<int, int>> lengths {};
        map<pair<int, int>, pair<int, int>> staple_counts {};
        set<pair<int, int>> staples {};
        map<pair<int, int>, pair<int, int>> transfs_counts {};
        set<pair<int, int>> transfs {};
        for (auto tracker: m_tracking) {
            auto info = tracker.first;
            auto counts = tracker.second;
            pair<int, int> length {info.num_linker_domains, info.num_central_domains};
            pair<int, int> num_staples {info.num_linker_staples, info.num_central_staples};
            pair<int, int> transf {info.disp_sum, info.rot_turns};
            lengths.insert(length);
            staples.insert(num_staples);
            transfs.insert(transf);
            if (lengths.find(length) == lengths.end()) {
                length_counts[length] = {counts.attempts, counts.accepts};
            }
            else {
                length_counts[length].first += counts.attempts;
                length_counts[length].second += counts.accepts;
            }
            if (staples.find(num_staples) == lengths.end()) {
                staple_counts[num_staples] = {counts.attempts, counts.accepts};
            }
            else {
                staple_counts[num_staples].first += counts.attempts;
                staple_counts[num_staples].second += counts.accepts;
            }
            if (transfs.find(transf) == transfs.end()) {
                transfs_counts[transf] = {counts.attempts, counts.accepts};
            }
            else {
                transfs_counts[transf].first += counts.attempts;
                transfs_counts[transf].second += counts.accepts;
            }
        }
        for (auto l: lengths) {
            *log_stream << "    Number of linker/central domains: " << l.first << "/" << l.second << "\n";
            int ats {length_counts[l].first};
            int acs {length_counts[l].second};
            double freq {static_cast<double>(acs) / ats};
            *log_stream << "        Attempts: " << ats << "\n";
            *log_stream << "        Accepts: " << acs << "\n";
            *log_stream << "        Frequency: " << freq << "\n";
        }
        *log_stream << "\n";
        for (auto st: staples) {
            *log_stream << "    Number of linker/central staples: " << st.first << "/" << st.second << "\n";
            int ats {staple_counts[st].first};
            int acs {staple_counts[st].second};
            double freq {static_cast<double>(acs) / ats};
            *log_stream << "        Attempts: " << ats << "\n";
            *log_stream << "        Accepts: " << acs << "\n";
            *log_stream << "        Frequency: " << freq << "\n";
        }
        *log_stream << "\n";
        for (auto t: transfs) {
            *log_stream << "    Sum/number of displacement/turns: " << t.first << "/" << t.second << "\n";
            int ats {transfs_counts[t].first};
            int acs {transfs_counts[t].second};
            double freq {static_cast<double>(acs) / ats};
            *log_stream << "        Attempts: " << ats << "\n";
            *log_stream << "        Accepts: " << acs << "\n";
            *log_stream << "        Frequency: " << freq << "\n";
        }
    }

    void CTCBLinkerRegrowthMCMovetype::reset_internal() {
        CTCBRegrowthMCMovetype::reset_internal();
        m_linker_endpoints.clear();
    }

    set<int> CTCBLinkerRegrowthMCMovetype::select_and_setup_segments(
            vector<Domain*>& linker1,
            vector<Domain*>& linker2,
            vector<Domain*>& central_segment) {

        // Pick region to be modified
        vector<Domain*> scaffold {m_origami_system.get_chain(
                m_origami_system.c_scaffold)};
        pair<int, int> endpoints {select_endpoints(scaffold.size(), 0, 3)};
        Domain* scaffold_start {scaffold[endpoints.first]};
        Domain* scaffold_end {scaffold[endpoints.second]};
        vector<Domain*> scaffold_domains {};
        if (m_origami_system.m_cyclic) {
            scaffold_domains = select_cyclic_segment(scaffold_start, scaffold_end);
        }
        else {
            scaffold_domains = select_linear_segment(scaffold_start, scaffold_end);
        }

        // Select endpoints (these will be inluded as central segment)
        // Exclude the terminal domains to ensure always linkers to regrow
        pair<int, int> internal_endpoints {select_internal_endpoints(
                scaffold_domains)};
        int start_di {internal_endpoints.first};
        int end_di {internal_endpoints.second};
        if (end_di < start_di) {
            std::swap(start_di, end_di);
        }
        
        // Get linker domains
        for (int i {0}; i != static_cast<int>(scaffold_domains.size()); i++) {
            if (i <= start_di) {
                linker1.push_back(scaffold_domains[i]);
            }
            if (i >= start_di and i <= end_di) {
                central_segment.push_back(scaffold_domains[i]);
            }
            if (i >= end_di) {
                linker2.push_back(scaffold_domains[i]);
            }
        }

        // Growing from external to central, order needs to reflect this
        std::reverse(linker1.begin(), linker1.end());

        return setup_fixed_end_biases(linker1, linker2, scaffold_domains);
    }

    set<int> CTCBLinkerRegrowthMCMovetype::setup_fixed_end_biases(
            vector<Domain*>& linker1,
            vector<Domain*>& linker2,
            vector<Domain*>& scaffold_domains) {

        // Add terminal constraint points
        int dir {scaffold_domains[1]->m_d - scaffold_domains[0]->m_d};
        Domain* linker1_endpoint {*scaffold_domains[0] + (-dir)};
        Domain* linker2_endpoint {*scaffold_domains.back() + dir};
        for (auto linker_endpoint: {linker1_endpoint, linker2_endpoint}) {
            m_linker_endpoints.push_back(linker_endpoint);
            if (linker_endpoint != nullptr) {
                m_constraintpoints.add_active_endpoint(linker_endpoint,
                        linker_endpoint->m_pos);
            }
        }

        // Find the rest
        vector<Domain*> linkers {};
        linkers.insert(linkers.end(), linker1.begin() + 1, linker1.end());
        linkers.insert(linkers.end(), linker2.begin() + 1, linker2.end());
        m_constraintpoints.calculate_constraintpoints(linkers, {});
        set<int> staples {m_constraintpoints.staples_to_be_regrown()};

        return staples;
    }

    pair<int, int> CTCBLinkerRegrowthMCMovetype::select_internal_endpoints(
            vector<Domain*> domains) {

        pair<int, int> internal_endpoints {select_endpoints(
                domains.size(), 1, 0)};

        return internal_endpoints;
    }

    bool CTCBLinkerRegrowthMCMovetype::domains_bound_externally(
            vector<Domain*> domains) {

        bool externally_bound {false};
        for (auto domain: domains) {
            if (domain->m_state != Occupancy::unbound) {
                Domain* bound_domain {domain->m_bound_domain};
                if (bound_domain->m_c == m_origami_system.c_scaffold) {
                    bool domain_in_range {find(domains.begin(), domains.end(),
                            bound_domain) != domains.end()};
                    if (not domain_in_range) {
                        externally_bound = true;
                        break;
                    }
                    else {
                        continue;
                    }
                }

                set<int> participating_chains {domain->m_c};
                if (not scan_for_external_scaffold_domain(bound_domain, domains,
                            participating_chains)) {
                    externally_bound = false;
                }
                else {
                    externally_bound = true;
                    break;
                }
            }
        }

        return externally_bound;
    }

    bool CTCBLinkerRegrowthMCMovetype::scan_for_external_scaffold_domain(
            Domain* domain,
            vector<Domain*> domains,
            set<int>& participating_chains) {

        bool externally_bound {false};
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
                    bool domain_in_range {find(domains.begin(), domains.end(),
                            bound_domain) != domains.end()};
                    if (not domain_in_range) {
                        externally_bound = true;
                        break;
                    }
                    else {
                        continue;
                    }
                }

                // Check if bound domain already in progress
                if (participating_chains.count(bound_domain->m_c) > 0) {
                    continue;
                }
                else {
                    externally_bound = scan_for_external_scaffold_domain(
                            bound_domain, domains, participating_chains);
                }
            }
        }

        return externally_bound;
    }

    void CTCBLinkerRegrowthMCMovetype::transform_segment(
            vector<Domain*> linker1,
            vector<Domain*> linker2,
            vector<Domain*> central_segment,
            vector<Domain*> central_domains) {

        // Try transformations until is not going to be immediately rejected
        bool regrowth_possible {false};
        while (not regrowth_possible) {

            // Translation component
            VectorThree disp {};
            for (int i {0}; i != 3; i++) {
                disp[i] = m_random_gens.uniform_int(0, m_max_disp);
            }

            // Rotation component
            // Select rotation center (from central scaffold domain positions)
            int center_di {m_random_gens.uniform_int(0, central_segment.size() - 1)};
            pair<int, int> center_key {central_segment[center_di]->m_c,
                    central_segment[center_di]->m_d};
            VectorThree center {m_prev_pos[center_key]};

            // Select axis and number of turns
            int axis_i {m_random_gens.uniform_int(0, 2)};
            VectorThree axis {utility::basis_vectors[axis_i]};
            int turns {m_random_gens.uniform_int(0, m_max_turns)};

            // Apply transformation
            bool transform_applied {apply_transformation(central_domains, disp,
                    center, axis, turns)};

            // Check if enough steps to reach endpoints
            if (transform_applied) {
                if (steps_less_than_distance(linker1, linker2)) {
                    reset_segment(central_domains, central_domains.size());
                }
                else {
                    regrowth_possible = true;
                    m_tracker.disp_sum = disp.sum();
                    m_tracker.rot_turns = turns;
                }
            }
        }
        write_config();
    }

    void CTCBLinkerRegrowthMCMovetype::reset_segment(vector<Domain*> segment,
            size_t last_di) {

        for (size_t di {0}; di != last_di; di++) {
            Domain* domain {segment[di]};
            m_origami_system.unassign_domain(*domain);
            m_assigned_domains.pop_back(); // HACK
        }
    }

    bool CTCBLinkerRegrowthMCMovetype::apply_transformation(
            vector<Domain*> central_domains,
            VectorThree disp,
            VectorThree center,
            VectorThree axis,
            int turns) {

        bool transform_applied {true};
        for (size_t di {0}; di != central_domains.size(); di++) {
            Domain* domain {central_domains[di]};
            pair<int, int> key {domain->m_c, domain->m_d};
            VectorThree pos {m_prev_pos[key]};
            VectorThree ore {m_prev_ore[key]};

            // Translation
            pos = pos + disp;

            // Rotation
            pos = pos.rotate(center, axis, turns);
            ore = ore.rotate({0, 0, 0}, axis, turns);

            // If position occupied by external domain, try again
            if (m_origami_system.position_occupancy(pos) == Occupancy::bound or
                    m_origami_system.position_occupancy(pos) == Occupancy::misbound) {
                reset_segment(central_domains, di);
                transform_applied = false;
                break;
            }
            else if (m_origami_system.position_occupancy(pos) == Occupancy::unbound) {
                Domain* unbound_domain {m_origami_system.unbound_domain_at(pos)};
                if (find(central_domains.begin(), central_domains.end(),
                            unbound_domain) == central_domains.end()) {
                    reset_segment(central_domains, di);
                    transform_applied = false;
                    break;
                }
            }
            m_origami_system.set_checked_domain_config(*domain, pos, ore);
            m_assigned_domains.push_back(key);
        }

        return transform_applied;
    }

    bool CTCBLinkerRegrowthMCMovetype::steps_less_than_distance(
            vector<Domain*> linker1,
            vector<Domain*> linker2) {

        bool steps_less {true};
        Domain* linker1_endpoint {m_linker_endpoints[0]};
        int dist1 {0};
        if (linker1_endpoint != nullptr) {
            dist1 = (linker1[0]->m_pos - linker1_endpoint->m_pos).abssum();
        }
        Domain* linker2_endpoint {m_linker_endpoints[1]};
        int dist2 {0};
        if (linker2_endpoint != nullptr) {
            dist2 = (linker2[0]->m_pos - linker2_endpoint->m_pos).abssum();
        }
        if (dist1 <= static_cast<int>(linker1.size()) and
            dist2 <= static_cast<int>(linker2.size())) {
            steps_less = false;
        }

        return steps_less;
    }

    void CTCBLinkerRegrowthMCMovetype::revert_transformation(
            vector<Domain*> central_domains) {

        for (auto domain: central_domains) {
            pair<int, int> key {domain->m_c, domain->m_d};
            VectorThree pos {m_old_pos[key]};
            VectorThree ore {m_old_ore[key]};
            m_origami_system.set_checked_domain_config(*domain, pos, ore);
            m_assigned_domains.push_back(key);
        }
        write_config();
    }

    pair<int, int> ClusteredCTCBLinkerRegrowth::select_internal_endpoints(
            vector<Domain*> domains) {

        int kernel {m_random_gens.uniform_int(1, domains.size() - 2)};
        int num_domains {static_cast<int>(domains.size())};
        bool segment_started {false};
        int start_di {1};
        int end_di {static_cast<int>(domains.size()) - 2};
        if (domains[kernel]->m_state == Occupancy::bound) {
            segment_started = true;
            for (int i {kernel - 1}; i != 0; i--) {
                if (domains[i]->m_state != Occupancy::bound) {
                    start_di = i + 1;
                    break;
                }
            }
        }
        for (int i {kernel + 1}; i != num_domains - 1; i++) {
            if (domains[i]->m_state == Occupancy::bound and not segment_started) {
                segment_started = true;
                start_di = i;
            }
            else if (domains[i]->m_state != Occupancy::bound and segment_started) {
                end_di = i - 1;
                break;
            }
        }
        pair<int, int> endpoints {start_di, end_di};
        if (not segment_started) {
            endpoints = CTCBLinkerRegrowthMCMovetype::select_internal_endpoints(
                    domains);
        }

        return endpoints;
    }
}
