// transform_movetypes.cpp

#include <map>

#include "utility.h"
#include "transform_movetypes.h"

namespace movetypes {

    using std::min;
    using std::map;

    using utility::Occupancy;

    LinkerRegrowthMCMovetype::LinkerRegrowthMCMovetype(
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
            int max_turns,
            unsigned int max_regrowth,
            unsigned int max_linker_length):
            MCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            CTRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_regrowth, max_regrowth),
            m_max_disp {max_disp},
            m_max_turns {max_turns},
            m_max_linker_length {max_linker_length} {
    }

    void LinkerRegrowthMCMovetype::write_log_summary(ostream* log_stream) {
        write_log_summary_header(log_stream);

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
            pair<int, int> length {info.num_linker_domains,
                    info.num_central_domains};
            pair<int, int> num_staples {info.num_linker_staples,
                    info.num_central_staples};
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
            *log_stream << "    Number of linker/central domains: " <<
                    l.first << "/" << l.second << "\n";
            int ats {length_counts[l].first};
            int acs {length_counts[l].second};
            double freq {static_cast<double>(acs) / ats};
            *log_stream << "        Attempts: " << ats << "\n";
            *log_stream << "        Accepts: " << acs << "\n";
            *log_stream << "        Frequency: " << freq << "\n";
        }
        *log_stream << "\n";
        for (auto st: staples) {
            *log_stream << "    Number of linker/central staples: " <<
                    st.first << "/" << st.second << "\n";
            int ats {staple_counts[st].first};
            int acs {staple_counts[st].second};
            double freq {static_cast<double>(acs) / ats};
            *log_stream << "        Attempts: " << ats << "\n";
            *log_stream << "        Accepts: " << acs << "\n";
            *log_stream << "        Frequency: " << freq << "\n";
        }
        *log_stream << "\n";
        for (auto t: transfs) {
            *log_stream << "    Sum/number of displacement/turns: " <<
                    t.first << "/" << t.second << "\n";
            int ats {transfs_counts[t].first};
            int acs {transfs_counts[t].second};
            double freq {static_cast<double>(acs) / ats};
            *log_stream << "        Attempts: " << ats << "\n";
            *log_stream << "        Accepts: " << acs << "\n";
            *log_stream << "        Frequency: " << freq << "\n";
        }
    }

    void LinkerRegrowthMCMovetype::add_tracker(bool accepted) {
        movetypes::add_tracker(m_tracker, m_tracking, accepted);
    }

    void LinkerRegrowthMCMovetype::reset_segment(vector<Domain*> segment,
            size_t last_di) {

        for (size_t di {0}; di != last_di; di++) {
            Domain* domain {segment[di]};
            m_origami_system.unassign_domain(*domain);
            m_assigned_domains.pop_back(); // HACK
        }
    }

    set<int> LinkerRegrowthMCMovetype::select_and_setup_segments(
            vector<Domain*>& linker1,
            vector<Domain*>& linker2,
            vector<Domain*>& central_segment) {

        // Reject moves that have central region bound externally
        bool externally_bound {true};

        // Could make max attempts settable
        int attempts {0};
        while (externally_bound and attempts != 10) {
            linker1.clear();
            linker2.clear();
            central_segment.clear();

            // Select region to be transformed
            central_segment = select_indices(m_scaffold, 1);

            // Select linker regions
            size_t linker1_length {static_cast<size_t>(
                        m_random_gens.uniform_int(2, m_max_linker_length + 1))};
            Domain* linker_d {central_segment.front()};
            while (linker1.size() != linker1_length and linker_d != nullptr and
                    linker_d != (*central_segment.back() + 2*m_dir)) {

                linker1.push_back(linker_d);
                linker_d = *linker_d + -m_dir;
            }
            if (m_origami_system.m_cyclic and linker1.size() == 1) {
                m_rejected = true;
                return {};
            }

            size_t linker2_length {static_cast<size_t>(
                        m_random_gens.uniform_int(2, m_max_linker_length + 1))};
            linker_d = central_segment.back();
            while (linker2.size() != linker2_length and linker_d != nullptr and
                    linker_d != (*linker1.back() + -m_dir)) {

                linker2.push_back(linker_d);
                linker_d = *linker_d + m_dir;
            }

            externally_bound = domains_bound_externally(central_segment);
            attempts++;
        }
        if (externally_bound) {
            m_rejected = true;
            return {};
        }

        return setup_fixed_end_biases(linker1, linker2);
    }

    set<int> LinkerRegrowthMCMovetype::setup_fixed_end_biases(
            vector<Domain*>& linker1,
            vector<Domain*>& linker2) {

        // Add terminal constraint points
        Domain* linker1_endpoint {*linker1.back() + -m_dir};
        Domain* linker2_endpoint {*linker2.back() + m_dir};
        int seg {0};
        for (auto linker_endpoint: {linker1_endpoint, linker2_endpoint}) {
            m_linker_endpoints.push_back(linker_endpoint);
            if (linker_endpoint != nullptr) {
                m_constraintpoints.add_active_endpoint(linker_endpoint,
                        linker_endpoint->m_pos, seg);
            }
            seg++;
        }

        // Find the rest
        vector<vector<Domain*>> linkers {{linker1.begin() + 1, linker1.end()},
                {linker2.begin() + 1, linker2.end()}};
        vector<int> dirs {-m_dir, m_dir};
        m_constraintpoints.calculate_constraintpoints(linkers, dirs, {});
        set<int> staples {m_constraintpoints.staples_to_be_regrown()};

        return staples;
    }

    bool LinkerRegrowthMCMovetype::domains_bound_externally(
            vector<Domain*> domains) {

        bool externally_bound {false};
        for (auto domain: domains) {
            if (domain->m_state != Occupancy::unbound) {
                Domain* bound_domain {domain->m_bound_domain};
                if (bound_domain->m_c == m_origami_system.c_scaffold) {
                    continue;
                }
                else {
                    set<int> participating_chains {domain->m_c};
                    if (scan_for_external_scaffold_domain(bound_domain, domains,
                                participating_chains)) {
                        externally_bound = true;
                        break;
                    }
                }
            }
        }

        return externally_bound;
    }

    bool LinkerRegrowthMCMovetype::scan_for_external_scaffold_domain(
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
                bool domain_on_scaffold {bound_domain->m_c ==
                        m_origami_system.c_scaffold};
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
                    if (externally_bound) {
                        break;
                    }
                }
            }
        }

        return externally_bound;
    }

    double LinkerRegrowthMCMovetype::transform_segment(
            vector<Domain*> linker1,
            vector<Domain*> linker2,
            vector<Domain*> central_segment,
            vector<Domain*> central_domains) {

        // Try transformations until is not going to be immediately rejected
        bool regrowth_possible {false};

        // Could make max attempts settable
        int attempts {0};
        double delta_e {0};
        // This probably doesn't strictly obey detailed balance
        while (not regrowth_possible and attempts != 1000) {

            // Translation component
            VectorThree disp {};
            for (int i {0}; i != 3; i++) {
                disp[i] = m_random_gens.uniform_int(-m_max_disp, m_max_disp);
            }

            // Rotation component
            // Select rotation center (from central scaffold domain positions)
            int center_di {m_random_gens.uniform_int(0,
                    central_segment.size() - 1)};
            pair<int, int> center_key {central_segment[center_di]->m_c,
                    central_segment[center_di]->m_d};
            VectorThree center {m_prev_pos[center_key]};

            // Select axis and number of turns
            int axis_i {m_random_gens.uniform_int(0, 2)};
            VectorThree axis {utility::basis_vectors[axis_i]};
            int turns {m_random_gens.uniform_int(0, m_max_turns)};

            // Apply transformation
            delta_e = apply_transformation(central_domains, disp,
                    center, axis, turns);

            // Check if enough steps to reach endpoints
            if (not m_transform_rejected) {
                if (steps_less_than_distance(linker1, linker2)) {
                    reset_segment(central_domains, central_domains.size());
                }
                else {
                    regrowth_possible = true;
                    m_tracker.disp_sum = disp.sum();
                    m_tracker.rot_turns = turns;
                }
            }
            attempts++;
        }
        if (not regrowth_possible) {
            m_rejected = true;
        }
        write_config();

        return delta_e;
    }

    double LinkerRegrowthMCMovetype::apply_transformation(
            vector<Domain*> central_domains,
            VectorThree disp,
            VectorThree center,
            VectorThree axis,
            int turns) {

        m_transform_rejected = false;
        double delta_e {0};
        for (size_t di {0}; di != central_domains.size(); di++) {
            Domain* domain {central_domains[di]};
            pair<int, int> key {domain->m_c, domain->m_d};
            VectorThree pos {m_prev_pos[key]};
            VectorThree ore {m_prev_ore[key]};

            // Rotation
            pos = pos.rotate(center, axis, turns);
            ore = ore.rotate({0, 0, 0}, axis, turns);

            // Translation
            pos = pos + disp;

            // If position occupied by external domain, try again
            if (m_origami_system.position_occupancy(pos) == Occupancy::bound or
                    m_origami_system.position_occupancy(pos) ==
                    Occupancy::misbound) {
                reset_segment(central_domains, di);
                m_transform_rejected = true;
                break;
            }
            else if (m_origami_system.position_occupancy(pos) ==
                    Occupancy::unbound) {
                Domain* unbound_domain {m_origami_system.unbound_domain_at(
                        pos)};
 
                bool scaffold_misbinding {domain->m_c ==
                        m_origami_system.c_scaffold and unbound_domain->m_c ==
                        m_origami_system.c_scaffold};
                bool new_binding_pair {find(central_domains.begin(),
                        central_domains.end(), unbound_domain) ==
                        central_domains.end()};
                if (not scaffold_misbinding and new_binding_pair) {
                    reset_segment(central_domains, di);
                    m_transform_rejected = true;
                    break;
                }
                else if (scaffold_misbinding and new_binding_pair) {
                    m_origami_system.check_domain_constraints(
                            *domain, pos, ore);
                    if (m_origami_system.m_constraints_violated) {
                        m_transform_rejected = true;
                        m_origami_system.m_constraints_violated = false;
                        reset_segment(central_domains, di);
                        break;
                    }
                }
            }
            delta_e += m_origami_system.set_checked_domain_config(*domain,
                    pos, ore);
            m_assigned_domains.push_back(key);
        }

        return delta_e;
    }

    bool LinkerRegrowthMCMovetype::steps_less_than_distance(
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

    ClusteredLinkerRegrowthMCMovetype::ClusteredLinkerRegrowthMCMovetype(
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
            int max_turns,
            unsigned int max_linker_length):
            MCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            CTRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, 0, 0),
            LinkerRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_disp, max_turns, 0,
                    max_linker_length) {
    }

    set<int> ClusteredLinkerRegrowthMCMovetype::select_and_setup_segments(
            vector<Domain*>& linker1,
            vector<Domain*>& linker2,
            vector<Domain*>& central_segment) {

        // Select a domain on the scaffold that is bound and then make central
        // region be all contiguously bound segments

        // Reject moves that have central region bound externally
        bool externally_bound {true};
        int attempts {0};
        while (externally_bound and attempts != 10) {
            if (m_origami_system.m_cyclic) {
                central_segment = cyclic_select_central_segment();
            }
            else {
                central_segment = linear_select_central_segment();
            }
            if (central_segment.size() == 0) {
                m_rejected = true;
                return {};
            }
            externally_bound = domains_bound_externally(central_segment);
            attempts++;
        }
        if (externally_bound) {
            m_rejected = true;
            return {};
        }

        // Get linker domains
        linker1.push_back(central_segment[0]);
        Domain* p_domain {*central_segment[0] + -1};
        if (p_domain != central_segment.back() and p_domain !=
                (*central_segment.back() + 1)) {
            while (p_domain != nullptr and
                    p_domain->m_state != Occupancy::bound and
                    linker1.size() != m_max_linker_length and
                    p_domain != (*central_segment.back() + 2)) {

                linker1.push_back(p_domain);
                p_domain = *p_domain + -1;
            }
        }
        if (linker1.size() == 1) {
            m_rejected = true;
            return {};
        }

        linker2.push_back(central_segment.back());
        Domain* n_domain {*central_segment.back() + 1};
        if (n_domain != linker1.back()) {
            while (n_domain != nullptr and
                    n_domain->m_state != Occupancy::bound and
                    linker2.size() != m_max_linker_length and
                    n_domain != (*linker1.back() + -1)) {

                linker2.push_back(n_domain);
                n_domain = *n_domain + 1;
            }
        }

        return setup_fixed_end_biases(linker1, linker2);
    }

    vector<Domain*> ClusteredLinkerRegrowthMCMovetype::linear_select_central_segment() {

        // Find a bound domain in the segment and set central cluster to be all
        // contiguous bound domains
        m_dir = 1;
        vector<Domain*> central_segment {};
        int kernel {m_random_gens.uniform_int(1, m_scaffold.size() - 2)};
        bool segment_started {false};
        Domain* kernel_domain {m_scaffold[kernel]};
        if (kernel_domain->m_state == Occupancy::bound) {
            segment_started = true;
            central_segment.push_back(kernel_domain);
            Domain* p_domain {*kernel_domain + -1};
            while (p_domain != nullptr and p_domain->m_state ==
                    Occupancy::bound and central_segment.size() !=
                    m_scaffold.size() - 1) {
                central_segment.push_back(p_domain);
                p_domain = *p_domain + -1;
            }
        }
        if (central_segment.size() == m_scaffold.size() - 1) {
            return {};
        }
        std::reverse(central_segment.begin(), central_segment.end());
        Domain* n_domain {*kernel_domain + 1};
        bool segment_ended {false};
        while (not segment_ended and n_domain != nullptr) {
            if (n_domain->m_state == Occupancy::bound and not
                    segment_started) {
                segment_started = true;
                central_segment.push_back(n_domain);
            }
            else if (n_domain->m_state == Occupancy::bound and
                    segment_started) {
                central_segment.push_back(n_domain);
            }
            else if (n_domain->m_state !=
                    Occupancy::bound and segment_started) {

                segment_ended = true;
            }
            n_domain = *n_domain + 1;
        }
        n_domain = *kernel_domain + -1;
        if (not segment_started) {
            while (not segment_ended and n_domain != nullptr) {
                if (n_domain->m_state == Occupancy::bound and not
                        segment_started) {
                    segment_started = true;
                    central_segment.push_back(n_domain);
                }
                else if (n_domain->m_state == Occupancy::bound and
                        segment_started) {
                    central_segment.push_back(n_domain);
                }
                else if (n_domain->m_state != Occupancy::bound and
                        segment_started) {

                    segment_ended = true;
                }
                n_domain = *n_domain + -1;
            }
            std::reverse(central_segment.begin(), central_segment.end());
        }

        return central_segment;
    }

    vector<Domain*> ClusteredLinkerRegrowthMCMovetype::cyclic_select_central_segment() {

        // Find a bound domain in the segment and set central cluster to be all
        // contiguous bound domains
        m_dir = 1;
        vector<Domain*> central_segment;
        int kernel {m_random_gens.uniform_int(1, m_scaffold.size() - 2)};
        bool segment_started {false};
        Domain* kernel_domain {m_scaffold[kernel]};
        if (kernel_domain->m_state == Occupancy::bound) {
            segment_started = true;
            central_segment.push_back(kernel_domain);
            Domain* p_domain {*kernel_domain + -1};
            while (p_domain->m_state == Occupancy::bound and
                    central_segment.size() != m_scaffold.size() - 1) {
                central_segment.push_back(p_domain);
                p_domain = *p_domain + -1;
            }
        }
        if (central_segment.size() == m_scaffold.size() - 1) {
            return {};
        }
        std::reverse(central_segment.begin(), central_segment.end());
        Domain* n_domain {*kernel_domain + 1};
        bool segment_ended {false};
        while (not segment_ended and n_domain != kernel_domain) {
            if (n_domain->m_state == Occupancy::bound and not
                    segment_started) {
                segment_started = true;
                central_segment.push_back(n_domain);
            }
            else if (n_domain->m_state == Occupancy::bound) {
                central_segment.push_back(n_domain);
            }
            else if (n_domain->m_state !=
                    Occupancy::bound and segment_started) {


                segment_ended = true;
            }
            n_domain = *n_domain + 1;
        }

        return central_segment;
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
            int max_turns,
            unsigned int max_regrowth,
            unsigned int max_linker_length):
            MCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            CTRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_regrowth, max_regrowth),
            LinkerRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_disp, max_turns,
                    max_regrowth, max_linker_length),
            RegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params),
            CBMCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            CTCBRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_regrowth) {
    }

    bool CTCBLinkerRegrowthMCMovetype::internal_attempt_move() {
        bool accepted {false};

        // Select region to be transformed and its linkers
        vector<Domain*> linker1 {};
        vector<Domain*> linker2 {};
        vector<Domain*> central_segment {};
        set<int> staples {select_and_setup_segments(linker1, linker2,
                central_segment)};
        if (m_rejected) {
            return accepted;
        }
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


        // Unassign domains
        vector<Domain*> linker1_to_unassign(linker1.begin() + 1, linker1.end());
        unassign_domains(linker1_to_unassign);
        vector<Domain*> linker2_to_unassign(linker2.begin() + 1, linker2.end());
        unassign_domains(linker2_to_unassign);
        for (auto c_i: staples) {
            unassign_domains(m_origami_system.get_chain(c_i));
        }

        // It would be nice if I didn't have to fully unassign them
        double delta_e {0};
        delta_e += unassign_domains(central_domains);

        delta_e += transform_segment(linker1, linker2, central_segment, central_domains);
        if (m_rejected) {
            return accepted;
        }
        m_bias *= exp(-delta_e);

        // Grow linkers
        for (auto linker: {linker1, linker2}) {
            grow_chain(linker);
            if (m_rejected) {
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

        // It would be nice if I didn't have to fully unassign them
        unassign_domains(central_domains);

        revert_transformation(central_domains);

        // Grow linkers
        for (auto linker: {linker1, linker2}) {
            grow_chain(linker);
        }

        // Reset modifier and test acceptance
        m_modifier = 1;
        accepted = test_cb_acceptance();

        return accepted;
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

    void CTCBLinkerRegrowthMCMovetype::reset_internal() {
        CTCBRegrowthMCMovetype::reset_internal();
        m_linker_endpoints.clear();
    }

    CTCBClusteredLinkerRegrowthMCMovetype::CTCBClusteredLinkerRegrowthMCMovetype(
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
            int max_turns,
            unsigned int max_linker_length):
            MCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            CTRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, 0, 0),
            LinkerRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_disp, max_turns, 0,
                    max_linker_length),
            RegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params),
            CBMCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            CTCBRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, 0),
            CTCBLinkerRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_disp, max_turns, 0,
                    max_linker_length),
            ClusteredLinkerRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_disp, max_turns,
                    max_linker_length) {
    }

    void CTCBClusteredLinkerRegrowthMCMovetype::reset_internal() {
        CTCBRegrowthMCMovetype::reset_internal();
        m_linker_endpoints.clear();
    }

    CTRGLinkerRegrowthMCMovetype::CTRGLinkerRegrowthMCMovetype(
            OrigamiSystem& origami_system,
            RandomGens& random_gens,
            IdealRandomWalks& ideal_random_walks,
            vector<OrigamiOutputFile*> config_files,
            string label,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params,
            int num_excluded_staples,
            int max_num_recoils,
            int max_c_attempts,
            int max_disp,
            int max_turns,
            unsigned int max_regrowth,
            unsigned int max_linker_length):
            MCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            CTRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_regrowth, max_regrowth),
            LinkerRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_disp, max_turns,
                    max_regrowth, max_linker_length),
            CTRGRegrowthMCMovetype(origami_system, random_gens,
                    ideal_random_walks, config_files, label, ops, biases,
                    params, num_excluded_staples, max_num_recoils,
                    max_c_attempts, max_regrowth) {
    }

    bool CTRGLinkerRegrowthMCMovetype::internal_attempt_move() {
        bool accepted {false};

        // Select region to be transformed and its linkers
        vector<Domain*> linker1 {};
        vector<Domain*> linker2 {};
        vector<Domain*> central_segment {};
        set<int> staples {select_and_setup_segments(linker1, linker2,
                central_segment)};
        if (m_rejected) {
            return accepted;
        }
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

        // Unassign domains
        m_regrow_ds = m_constraintpoints.domains_to_be_regrown();
        m_regrow_ds.insert(m_regrow_ds.begin(), linker1.front());
        m_delta_e += unassign_and_save_domains();
        m_delta_e += unassign_and_save_domains(central_domains);

        // Transform central segment
        m_delta_e += transform_segment(linker1, linker2, central_segment, central_domains);
        if (m_rejected) {
            return accepted;
        }

        // Grow linkers
        m_delta_e += recoil_regrow();
        if (m_rejected) {
            accepted = false;
            return accepted;
        }
        add_external_bias();

        // Calculate new weights
        setup_for_calc_new_weights();
        unassign_and_save_domains();
        m_constraintpoints.reset_active_endpoints();
        calc_weights();

        // Calculate old weights
        unassign_domains();
        unassign_and_save_domains(central_domains);
        revert_transformation(central_domains);
        m_constraintpoints.reset_active_endpoints();
        calc_old_c_opens();
        setup_for_calc_old_weights();
        unassign_and_save_domains();
        for (auto d: central_domains) {
            pair<int, int> key {d->m_c, d->m_d};
            m_modified_domains.push_back(key);
        }
        m_constraintpoints.reset_active_endpoints();
        calc_weights();

        // Test acceptance
        accepted = test_rg_acceptance();

        return accepted;
    }

    void CTRGLinkerRegrowthMCMovetype::revert_transformation(
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

    void CTRGLinkerRegrowthMCMovetype::reset_internal() {
        CTRGRegrowthMCMovetype::reset_internal();
        m_linker_endpoints.clear();
    }
}
