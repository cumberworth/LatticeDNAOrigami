// cb_movetypes.cpp

#include "algorithm"
#include "utility.h"
#include "movetypes.h"
#include "cb_movetypes.h"

namespace CBMovetypes {

    using namespace Movetypes;
    using namespace CBMovetypes;
    using namespace Utility;

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
                    //bool domains_complementary {m_origami_system.check_domains_complementary(
                    //        domain, *unbound_domain)};
                    VectorThree o_new;
                    //double bfactor;
                    //if (domains_complementary) {
                    //    o_new = -unbound_domain->m_ore;
                    //    bfactor = 1;
                    //}
                    //else {
                    //    o_new = {0, 0, 0};
                    //    bfactor = 6;
                    //}
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

    void CBMCMovetype::select_and_set_config(vector<Domain*> domains, int i) {
        Domain* domain {domains[i]};
        Domain* prev_domain {domains[i - 1]};
        VectorThree p_prev {prev_domain->m_pos};

        // Calculate weights
        vector<pair<VectorThree, VectorThree>> configs {};
        vector<double> bfactors {};
        calc_biases(*domain, p_prev, configs, bfactors);
        vector<double> weights {calc_bias(bfactors, domain, configs, p_prev, domains)};
        if (m_rejected) {
            return;
        }

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

    void CBMCMovetype::select_and_set_new_config(
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

    //double CBMCMovetype::set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old) {
    //    bool domains_complementary {m_origami_system.check_domains_complementary(
    //            growth_domain_new, growth_domain_old)};
    //    VectorThree o_new;
    //    double bias {1};
    //    if (domains_complementary) {
    //        o_new = -growth_domain_old.m_ore;
    //    }
    //    else {
    //        o_new = select_random_orientation();
    //        bias *= 6;
    //    }
    //    double delta_e {0};
    //    delta_e += m_origami_system.set_domain_config(growth_domain_new,
    //            growth_domain_old.m_pos, o_new);
    //    if (m_origami_system.m_constraints_violated) {
    //        m_rejected = true;
    //    }
    //    else {
    //        bias *= exp(-delta_e);
    //        m_bias *= bias;
    //        pair<int, int> key {growth_domain_new.m_c, growth_domain_new.m_d};
    //        m_assigned_domains.push_back(key);
    //    }
    //
    //    return delta_e;
    //}

    double CBMCMovetype::set_old_growth_point(Domain& growth_domain_new, Domain& growth_domain_old) {
        pair<int, int> key {growth_domain_new.m_c, growth_domain_new.m_d};
        VectorThree o_old {m_old_ore[key]};
        double delta_e {0};
        delta_e += m_origami_system.set_checked_domain_config(growth_domain_new,
                growth_domain_old.m_pos, o_old);
        //bool domains_complementary {m_origami_system.check_domains_complementary(
        //        growth_domain_new, growth_domain_old)};
        //double bias {1};
        //if (not domains_complementary) {
        //    bias *= 6;
        //}
        //bias *= exp(-delta_e);
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

    vector<pair<Domain*, Domain*>> CBMCMovetype::find_bound_domains(
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
            cout << "d\n";
            throw OrigamiMisuse {};
        }

        return bound_domains;
    }

    pair<Domain*, Domain*> CBMCMovetype::select_old_growthpoint(
            vector<pair<Domain*, Domain*>> bound_domains) {

        int bound_domain_index {m_random_gens.uniform_int(0, bound_domains.size() - 1)};
        Domain* growth_domain_new {bound_domains[bound_domain_index].first};
        Domain* growth_domain_old {bound_domains[bound_domain_index].second};
        return {growth_domain_new, growth_domain_old};
    }

    void CBMCMovetype::update_bias(int sign) {
        double total_bias {m_system_bias.calc_bias()};
        m_bias *= exp(-sign*total_bias);
    }

    bool CBStapleRegrowthMCMovetype::attempt_move() {
        bool accepted;

        // No staples to regrow
        if (m_origami_system.num_staples() == 0) {
            accepted = false;
            return accepted;
        }

        // Select a staple to regrow
        int c_i_index {m_random_gens.uniform_int(1, m_origami_system.num_staples())};
        vector<Domain*> selected_chain {m_origami_system.get_chains()[c_i_index]};

        // Reject if staple is connector
        if (staple_is_connector(selected_chain)) {
            accepted = false;
            return accepted;
        }

        // Select growth points on chains
        pair<Domain*, Domain*> growthpoint {select_new_growthpoint(selected_chain)};

        auto bound_domains = find_bound_domains(selected_chain);
        unassign_domains(selected_chain);

        // Grow staple
        set_growthpoint_and_grow_staple(growthpoint, selected_chain);
        if (m_rejected) {
            accepted = false;
            return accepted;
        }

        //HACK
        update_bias(+1);
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
        //HACK
        update_bias(+1);
        accepted = test_cb_acceptance();
        return accepted;
    }

    void CBStapleRegrowthMCMovetype::set_growthpoint_and_grow_staple(
            pair<Domain*, Domain*> growthpoint, vector<Domain*> selected_chain) {

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
            select_and_set_config(domains, i);
            if (m_rejected) {
                break;
            }
        }
    }

    vector<double> CBStapleRegrowthMCMovetype::calc_bias(vector<double> bfactors,
            Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree,
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

    bool CTCBScaffoldRegrowthMCMovetype::attempt_move() {
        bool accepted;

        m_regrow_old = false;
        vector<Domain*> scaffold {m_origami_system.get_chain(
                m_origami_system.c_scaffold)};
        vector<Domain*> scaffold_domains {select_indices(scaffold)};

        m_constraintpoints.calculate_constraintpoints(scaffold_domains);
        set<int> staples {m_constraintpoints.staples_to_be_regrown()};
        
        // Unassign staples except those directly and indirectly bound to external scaffold domains
        vector<Domain*> scaffold_domains_to_unassign(scaffold_domains.begin() + 1,
                scaffold_domains.end());
        unassign_domains(scaffold_domains_to_unassign);
        for (auto c_i: staples) {
            unassign_domains(m_origami_system.get_chain(c_i));
        }

        // Grow scaffold and staples
        if (m_constraintpoints.is_growthpoint(scaffold_domains[0])) {
            grow_staple_and_update_endpoints(scaffold_domains[0]);
            if (m_rejected) {
                accepted = false;
                return accepted;
            }
        }
        grow_chain(scaffold_domains);
        if (m_rejected) {
            accepted = false;
            return accepted;
        }

        //HACK
        update_bias(+1);
        // Regrow in old conformation
        setup_for_regrow_old();
        m_constraintpoints.reset_active_endpoints();

        // Unassign staples except those directly and indirectly bound to external scaffold domains
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
        //HACK
        update_bias(+1);
        accepted = test_cb_acceptance();
        return accepted;
    }

    void CTCBScaffoldRegrowthMCMovetype::reset_internal() {
        CBMCMovetype::reset_internal();
        m_constraintpoints.reset_internal();
    }

    vector<Domain*> CTCBScaffoldRegrowthMCMovetype::select_indices(
            vector<Domain*> segment) {

        pair<Domain*, Domain*> endpoints {select_endpoints(segment, 1)};
        vector<Domain*> domains {};
        if (m_origami_system.m_cyclic) {
            domains = select_cyclic_segment(endpoints.first, endpoints.second);
        }
        else {
            domains = select_linear_segment(endpoints.first, endpoints.second);
        }

        // If end domain is end of chain, no endpoint
        Domain* endpoint_domain {(*endpoints.second) + m_dir};
        if (endpoint_domain != nullptr) {
            m_constraintpoints.add_active_endpoint(endpoint_domain, endpoint_domain->m_pos);
        }

        return domains;
    }

    pair<Domain*, Domain*> CTCBScaffoldRegrowthMCMovetype::select_endpoints(
            vector<Domain*> domains, int min_size) {

        Domain* start_domain {domains[m_random_gens.uniform_int(0, domains.size() - 1)]};
        Domain* end_domain {domains[m_random_gens.uniform_int(0, domains.size() - 1)]};
        while (std::abs(start_domain->m_d - end_domain->m_d) < min_size) {
            end_domain = domains[m_random_gens.uniform_int(0, domains.size() - 1)];
        }

        return {start_domain, end_domain};
    }

    vector<Domain*> CTCBScaffoldRegrowthMCMovetype::select_cyclic_segment(
            Domain* start_domain, Domain* end_domain) {

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

    vector<Domain*> CTCBScaffoldRegrowthMCMovetype::select_linear_segment(
            Domain* start_domain, Domain* end_domain) {

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

    void CTCBScaffoldRegrowthMCMovetype::grow_chain(vector<Domain*> domains) {
        for (size_t i {1}; i != domains.size(); i++) {
            select_and_set_config(domains, i);
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

    void CTCBScaffoldRegrowthMCMovetype::grow_staple_and_update_endpoints(
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

    vector<double> CTCBScaffoldRegrowthMCMovetype::calc_bias(
            vector<double> bfactors,
            Domain* domain,
            vector<pair<VectorThree, VectorThree>>& configs,
            VectorThree p_prev,
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
                if (not (binding_same_chain or endpoint)) {
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

    bool CTCBLinkerRegrowthMCMovetype::attempt_move() {
        bool accepted;
        m_regrow_old = false;

        // Pick region to be modified
        vector<Domain*> scaffold {m_origami_system.get_chain(
                m_origami_system.c_scaffold)};
        pair<Domain*, Domain*> endpoints {select_endpoints(scaffold, 3)};
        vector<Domain*> scaffold_domains {};
        if (m_origami_system.m_cyclic) {
            scaffold_domains = select_cyclic_segment(endpoints.first, endpoints.second);
        }
        else {
            scaffold_domains = select_linear_segment(endpoints.first, endpoints.second);
        }

        // Pick linker regions/central region
        vector<Domain*> linker1 {};
        vector<Domain*> linker2 {};
        vector<Domain*> central_segment {};
        select_linkers(scaffold_domains, linker1, linker2, central_segment);

        // Reject moves that have central region bound externally
        if (domains_bound_externally(central_segment)) {
            accepted = false;
            return accepted;
        }

        // Calculate topology constraints
        vector<Domain*> linkers {linker1};
        linkers.insert(linkers.end(), linker2.begin(), linker2.end());
        m_constraintpoints.calculate_constraintpoints(linkers); // Deal with terminal domains properly
        set<int> staples {m_constraintpoints.staples_to_be_regrown()};
        
        // Unassign domains
        vector<Domain*> linker1_domains_to_unassign(linker1.begin() + 1,
                linker1.end());
        unassign_domains(linker1_domains_to_unassign);
        vector<Domain*> linker2_domains_to_unassign(linker2.begin() + 1,
                linker2.end());
        unassign_domains(linker2_domains_to_unassign);
        for (auto c_i: staples) {
            unassign_domains(m_origami_system.get_chain(c_i));
        }
        set<int> central_staples {find_staples(central_segment)};
        vector<Domain*> central_domains {central_segment};
        for (auto staple: central_staples) {
            vector<Domain*> staple_domains {m_origami_system.get_chain(staple)};
            for (auto domain: staple_domains) {
                central_domains.push_back(domain);
            }
        }

        // Consider a more efficient method for this
        unassign_domains(central_domains); 

        // Select and apply transformation for central segment
        bool regrowth_possible {false};
        while (not regrowth_possible) {

            // Translation component
            VectorThree disp {};
            for (int i {0}; i != 3; i++) {
                disp[i] = m_random_gens.uniform_int(0, m_params.m_max_displacement);
            }

            // Rotation component
            // Select rotation center (from central scaffold domain positions)
            int center_di {m_random_gens.uniform_int(0, central_segment.size() - 1)};
            pair<int, int> center_key {central_segment[center_di]->m_c,
                    central_segment[center_di]->m_d};
            VectorThree center {m_prev_pos[center_key]};

            // Select axis and number of turns
            int axis_i {m_random_gens.uniform_int(0, 2)};
            VectorThree axis {basis_vectors[axis_i]};
            int turns {m_random_gens.uniform_int(0, 3)};

            // Apply transformation
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

            // If regrowth not possible, reset central domain
            if (transform_applied) {
                int dist1 {(linker1[0]->m_pos - central_segment[0]->m_pos).abssum()};
                int dist2 {(linker2[0]->m_pos - central_segment.back()->m_pos).abssum()};
                if (dist1 <= static_cast<int>(linker1.size()) and
                        dist2 <= static_cast<int>(linker2.size())) {
                    regrowth_possible = true;
                }
                else {
                    reset_segment(central_domains, central_domains.size());
                }
            }
        }

        // Add terminal constraint points (consider only adding once if central seg is only 1 domain)
        m_constraintpoints.add_active_endpoint(central_segment[0],
                central_segment[0]->m_pos);
        if (central_segment[0]->m_d != central_segment.back()->m_d) {
            m_constraintpoints.add_active_endpoint(central_segment.back(),
                    central_segment.back()->m_pos);
        }

        // Grow linkers
        if (m_constraintpoints.is_growthpoint(linker1[0])) {
            grow_staple_and_update_endpoints(linker1[0]);
            if (m_rejected) {
                accepted = false;
                return accepted;
            }
        }
        grow_chain(linker1);
        if (m_rejected) {
            accepted = false;
            return accepted;
        }
        if (m_constraintpoints.is_growthpoint(linker2[0])) {
            grow_staple_and_update_endpoints(linker2[0]);
            if (m_rejected) {
                accepted = false;
                return accepted;
            }
        }
        grow_chain(linker2);
        if (m_rejected) {
            accepted = false;
            return accepted;
        }
        m_origami_system.check_all_constraints();

        //HACK
        update_bias(+1);
        // Regrow in old conformation
        setup_for_regrow_old();
        m_constraintpoints.reset_active_endpoints();

        // Unassign staples except those directly and indirectly bound to external scaffold domains
        unassign_domains(linker1_domains_to_unassign);
        unassign_domains(linker2_domains_to_unassign);
        for (auto c_i: staples) {
            unassign_domains(m_origami_system.get_chain(c_i));
        }

        // Consider a more efficient method for this
        unassign_domains(central_domains); 

        // Move central region back to original position
        for (auto domain: central_domains) {
            pair<int, int> key {domain->m_c, domain->m_d};
            VectorThree pos {m_old_pos[key]};
            VectorThree ore {m_old_ore[key]};
            m_origami_system.set_checked_domain_config(*domain, pos, ore);
            m_assigned_domains.push_back(key);
        }

        // Grow linkers
        if (m_constraintpoints.is_growthpoint(linker1[0])) {
            grow_staple_and_update_endpoints(linker1[0]);
        }
        grow_chain(linker1);
        if (m_constraintpoints.is_growthpoint(linker2[0])) {
            grow_staple_and_update_endpoints(linker2[0]);
        }
        grow_chain(linker2);
        m_origami_system.check_all_constraints();

        // Reset modifier and test acceptance
        m_modifier = 1;
        //HACK
        update_bias(+1);
        accepted = test_cb_acceptance();
        m_origami_system.check_all_constraints();
        return accepted;
    }

    void CTCBLinkerRegrowthMCMovetype::reset_segment(vector<Domain*> segment,
            size_t last_di) {

        for (size_t di {0}; di != last_di; di++) {
            Domain* domain {segment[di]};
            m_origami_system.unassign_domain(*domain);
            m_assigned_domains.pop_back(); // HACK
        }
    }

    void CTCBLinkerRegrowthMCMovetype::select_linkers(vector<Domain*> domains,
            vector<Domain*>& linker1, vector<Domain*>& linker2,
            vector<Domain*>& central_segment) {

        // Select endpoints (these will be inluded as central segment)
        // Exclude the terminal domains to ensure always linkers to regrow
        int start_di {m_random_gens.uniform_int(1, domains.size() - 2)};
        int end_di {m_random_gens.uniform_int(1, domains.size() - 2)};

        if (end_di < start_di) {
            std::swap(start_di, end_di);
        }
        
        for (int i {0}; i != static_cast<int>(domains.size()); i++) {
            if (i < start_di) {
                linker1.push_back(domains[i]);
            }
            else if (i >= start_di and i <= end_di) {
                central_segment.push_back(domains[i]);
            }
            else {
                linker2.push_back(domains[i]);
            }
        }

        // Growing from external to central, order needs to reflect this
        std::reverse(linker2.begin(), linker2.end());
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
}
