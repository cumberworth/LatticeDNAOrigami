// movetypes.cpp

#include <algorithm>
#include <random>
#include <set>
#include <cmath>

#include "movetypes.h"

namespace movetypes {

    using std::fmin;
    using std::min;
    using std::set;
    using std::find;
    using utility::Occupancy;
    using utility::OrigamiMisuse;

	MCMovetype::MCMovetype(
                OrigamiSystem& origami_system,
                RandomGens& random_gens,
                IdealRandomWalks& ideal_random_walks,
                vector<OrigamiOutputFile*> config_files,
                string label,
                SystemOrderParams& ops,
                SystemBiases& biases,
                InputParameters& params) :

        m_origami_system {origami_system},
        m_random_gens {random_gens},
        m_ideal_random_walks {ideal_random_walks},
        m_config_files {config_files},
        m_label {label},
        m_ops {ops},
        m_biases {biases},
        m_params {params},
        m_config_output_freq {params.m_vtf_output_freq},
        m_max_total_staples {params.m_max_total_staples},
        m_max_type_staples {params.m_max_type_staples} {
	}

    bool MCMovetype::attempt_move(long long int step) {
        reset_internal();
        m_step = step;
        write_config();
        m_general_tracker.attempts++;
        bool accepted {internal_attempt_move()};
        m_general_tracker.accepts += accepted;
        add_tracker(accepted);

        return accepted;
    }

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

    void MCMovetype::write_log_summary_header(ostream* log_stream) {
        *log_stream << "Movetype: " << get_label() << "\n";
        int attempts {get_attempts()};
        int accepts {get_accepts()};
        double freq {static_cast<double>(accepts) / attempts};
        *log_stream << "    Attempts: " << attempts << "\n";
        *log_stream << "    Accepts: " << accepts << "\n";
        *log_stream << "    Frequency: " << freq << "\n";
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
        VectorThree vec {utility::vectors[m_random_gens.uniform_int(0, 5)]};
        return p_prev + vec;
    }

    VectorThree MCMovetype::select_random_orientation() {
        VectorThree vec {utility::vectors[m_random_gens.uniform_int(0, 5)]};
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

    set<int> MCMovetype::find_staples(vector<Domain*> domains) {
        set<int> staples {};
        for (auto domain: domains) {
            Domain* bound_domain {domain->m_bound_domain};
            if (bound_domain != nullptr and bound_domain->m_c != 0) {
                scan_for_scaffold_domain(bound_domain, staples);
            }
        }

        return staples;
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

    void MCMovetype::write_config() {
        if (m_config_output_freq != 0 and m_step % m_config_output_freq == 0) {
            for (auto file: m_config_files) {
                file->write(0);
            }
        }
    }

    void MCMovetype::reset_internal() {
        m_modified_domains.clear();
        m_assigned_domains.clear();
        m_added_chains.clear();
        m_prev_pos.clear();
        m_prev_ore.clear();
        m_rejected = false;
        m_modifier = 1;
    }

    string MCMovetype::get_label() {
        return m_label;
    }

    int MCMovetype::get_attempts() {
        return m_general_tracker.attempts;
    }

    int MCMovetype::get_accepts() {
        return m_general_tracker.accepts;
    }

    vector<domainPairT> MCMovetype::find_bound_domains(
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

	IdentityMCMovetype::IdentityMCMovetype(
                OrigamiSystem& origami_system,
                RandomGens& random_gens,
                IdealRandomWalks& ideal_random_walks,
                vector<OrigamiOutputFile*> config_files,
                string label,
                SystemOrderParams& ops,
                SystemBiases& biases,
                InputParameters& params) :
        MCMovetype(origami_system, random_gens, ideal_random_walks,
                config_files, label, ops, biases, params) {
	}

	RegrowthMCMovetype::RegrowthMCMovetype(
                OrigamiSystem& origami_system,
                RandomGens& random_gens,
                IdealRandomWalks& ideal_random_walks,
                vector<OrigamiOutputFile*> config_files,
                string label,
                SystemOrderParams& ops,
                SystemBiases& biases,
                InputParameters& params) :
        MCMovetype(origami_system, random_gens, ideal_random_walks,
                config_files, label, ops, biases, params) {
	}

    double RegrowthMCMovetype::set_growth_point(
            Domain& growth_domain_new,
            Domain& growth_domain_old) {

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
            write_config();
        }

        return delta_e;
    }

    void RegrowthMCMovetype::grow_staple(int d_i_index, vector<Domain*> selected_chain) {
        // Grow staple in both directions out from growth point
        // The indices passed to grow chain should include the growth point

        // Grow in three prime direction (staple domains increase in 3' direction)
        auto first_iter3 {selected_chain.begin() + d_i_index};
        auto last_iter3 {selected_chain.end()};
        vector<Domain*> domains_three_prime {first_iter3, last_iter3};
        if (domains_three_prime.size() > 1) {
            grow_chain(domains_three_prime);
        }
        if (m_rejected) {
            return;
        }

        // Grow in five prime direction
        auto first_iter5 {selected_chain.begin()};
        auto last_iter5 {selected_chain.begin() + d_i_index + 1};
        vector<Domain*> domains_five_prime {first_iter5, last_iter5};
        std::reverse(domains_five_prime.begin(), domains_five_prime.end());
        if (domains_five_prime.size() > 1) {
            grow_chain(domains_five_prime);
        }
        if (m_rejected) {
            return;
        }
    }

    pair<Domain*, Domain*> RegrowthMCMovetype::select_new_growthpoint(
            vector<Domain*> selected_chain) {

        int growth_di_new {m_random_gens.uniform_int(0, selected_chain.size() -
                1)};
        Domain* growth_d_new {selected_chain[growth_di_new]};
        Domain* growth_d_old {select_random_domain()};
        while (growth_d_old->m_c == growth_d_new->m_c) {
            growth_d_old = select_random_domain();
        }

        return {growth_d_new, growth_d_old};
    }

    domainPairT RegrowthMCMovetype::select_old_growthpoint(
            vector<domainPairT> bound_domains) {

        int bound_domain_index {m_random_gens.uniform_int(0,
                bound_domains.size() - 1)};
        Domain* growth_domain_new {bound_domains[bound_domain_index].first};
        Domain* growth_domain_old {bound_domains[bound_domain_index].second};
        return {growth_domain_new, growth_domain_old};
    }

    int RegrowthMCMovetype::num_bound_staple_domains(vector<Domain*> staple) {
        int num_bd {0};
        for (auto d: staple) {
            if (d->m_state == Occupancy::bound or
                    d->m_state == Occupancy::misbound) {
                num_bd++;
            }
        }

        return num_bd;
    }

    CTRegrowthMCMovetype::CTRegrowthMCMovetype(
                OrigamiSystem& origami_system,
                RandomGens& random_gens,
                IdealRandomWalks& ideal_random_walks,
                vector<OrigamiOutputFile*> config_files,
                string label,
                SystemOrderParams& ops,
                SystemBiases& biases,
                InputParameters& params,
                int num_excluded_staples,
                int max_regrowth) :
        MCMovetype(origami_system, random_gens, ideal_random_walks,
                config_files, label, ops, biases, params),
        m_num_excluded_staples {num_excluded_staples},
        m_max_regrowth {static_cast<unsigned int>(max_regrowth)} {

            m_scaffold = m_origami_system.get_chain(m_origami_system.c_scaffold);
    }

    void CTRegrowthMCMovetype::reset_internal() {
        m_constraintpoints.reset_internal();
        m_excluded_staples.clear();
    }

    void CTRegrowthMCMovetype::sel_excluded_staples() {
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
    }

    vector<Domain*> CTRegrowthMCMovetype::select_indices(
            vector<Domain*> segment, unsigned int min_length, int seg) {

        unsigned int seg_length {static_cast<unsigned int>(segment.size())};
        unsigned int max_length {min(seg_length, m_max_regrowth)};
        unsigned int sel_length {static_cast<unsigned int>(
                m_random_gens.uniform_int(min_length, max_length))};
        int start_i {m_random_gens.uniform_int(0, seg_length - 1)};

        // Select direction of regrowth
        m_dir = m_random_gens.uniform_int(0, 1);
        if (m_dir == 0) {
            m_dir = -1;
        }

        // Add domains until length reached or end of segment reached
        Domain* cur_domain {segment[start_i]};
        vector<Domain*> domains {};
        while (cur_domain != nullptr and domains.size() != sel_length) {
            domains.push_back(cur_domain);
            cur_domain = (*cur_domain) + m_dir;
        }
        if (domains.size() < min_length) {
            domains = select_indices(segment, min_length, seg);
            return domains;
        }

        // If end domain is end of chain, no endpoint
        if (cur_domain != nullptr) {
            m_constraintpoints.add_active_endpoint(cur_domain,
                    cur_domain->m_pos, seg);
        }

        return domains;
    }

    bool CTRegrowthMCMovetype::excluded_staples_bound() {
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

        return bound_to_system;
    }

}
