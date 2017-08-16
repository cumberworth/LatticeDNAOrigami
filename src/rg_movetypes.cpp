// rg_movetypes.cpp

#include <rg_movetypes.h>

namespace movetypes {

    using utility::all_pairs;

    CTScaffoldRG::CTScaffoldRG(
            OrigamiSystem& origami_system,
            RandomGens& random_gens,
            IdealRandomWalks& ideal_random_walks,
            vector<OrigamiOutputFile*> config_files,
            string label,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params,
            int num_excluded_staples,
            int max_recoils):

        MCMovetype(origami_system, random_gens, ideal_random_walks,
                config_files, label, ops, biases, params),
        CTRegrowthMCMovetype(num_excluded_staples),
        m_all_configs {all_pairs<VectorThree>(utility::vectors)},
        m_all_cis (m_all_configs.size()),
        m_max_recoils {max_recoils} {

        std::iota(m_all_cis.begin(), m_all_cis.end(), 0);
    }

    void CTScaffoldRG::write_log_summary(ostream*) {
    }

    void CTScaffoldRG::reset_internal() {
        MCMovetype::reset_internal();
        m_sel_scaf_doms.clear();
        m_regrow_ds.clear();
        m_c_attempts_q.clear();
    }

    bool CTScaffoldRG::internal_attempt_move() {
        bool accepted {false};

        // Select scaffold indices and excluded staples
        m_sel_scaf_doms = select_indices(m_scaffold);
        sel_excluded_staples();

        // Setup up constraint points
        m_constraintpoints.calculate_constraintpoints(m_sel_scaf_doms,
                m_excluded_staples);

        // If scaffold is cyclic and the whole range is selected, the endpoints
        // on the first domain are needed, otherwise remove
        if (not m_origami_system.m_cyclic and m_sel_scaf_doms.size() !=
                m_origami_system.get_chain(0).size()) {

            m_constraintpoints.remove_active_endpoint(m_sel_scaf_doms[0]);
        }

        // Get domains to be regrown
        m_regrow_ds = m_constraintpoints.domains_to_be_regrown();
        //m_tracker.num_staples = staples.size();

        // Recoil regrow
        unassign_domains();
        recoil_regrow();

        // Reject if constraints disobeyed
        if (m_rejected) {
            accepted = false;
            return accepted;
        }

        // Reject if excluded staples not bound to system
        m_rejected = excluded_staples_bound();
        if (m_rejected) {
            accepted = false;
            return accepted;
        }

        // Calculate new weights
        add_external_bias();
        unassign_domains();
        calc_rg_weights();

        // Calculate old weights
        setup_for_regrow_old();
        unassign_domains();
        calc_rg_weights();
        add_external_bias();

        // Test acceptance
        accepted = test_rg_acceptance();

        return accepted;
    }

    void CTScaffoldRG::add_tracker(bool) {
    }

    void CTScaffoldRG::unassign_domains() {
        for (auto domain: m_regrow_ds) {
            pair<int, int> key {domain->m_c, domain->m_d};
            m_prev_pos[key] = domain->m_pos;
            m_prev_ore[key] = domain->m_ore;
            m_modified_domains.push_back(key);

            // Would be double counting to include delta_e in bias
            m_origami_system.unassign_domain(*domain);
        }
        write_config();
        update_external_bias(); // TODO you sure you want this here?
    }

    void CTScaffoldRG::update_external_bias() {
        m_ops.update_move_params();
        m_biases.calc_move();
    }

    void CTScaffoldRG::recoil_regrow() {
        bool growth_q_empty {false};
        bool recoil_q_empty {false};
        Domain* cur_d {m_regrow_ds.front()};
        m_regrow_ds.pop_front();
        deque<Domain*> recoil_ds {};
        int recoils {0};
        deque<int> c_attempts_q {};
        int c_attempts {0};
        int cur_c_max_attempts;
        bool growthpoint {m_constraintpoints.is_growthpoint(cur_d)};
        vector<int> cur_avail_cis;
        VectorThree prev_p {m_sel_scaf_doms[0]->m_pos};
        int dir {cur_d->m_d - m_sel_scaf_doms[0]->m_d};
        if (growthpoint) {
            cur_c_max_attempts = 1;
            cur_avail_cis = {};
            dir = m_regrow_ds.front()->m_d - cur_d->m_d;
        }
        else {
            cur_c_max_attempts = m_max_c_attempts;
            cur_avail_cis = m_all_cis;
            prev_p = (*cur_d + dir)->m_pos;
        }
        while (not (growth_q_empty or recoil_q_empty)) {
            bool c_open {false};
            configT config;
            while (not c_open and c_attempts != m_max_c_attempts) {
                c_attempts++;
                config = select_trial_config(cur_d, cur_avail_cis, growthpoint,
                        prev_p);
                c_open = test_config_open(config); 
            }
            if (c_open) {
                m_origami_system.set_domain_config(*cur_d, config.first,
                        config.second);
                if (recoils == m_max_recoils) {
                    Domain* confirmed_d {recoil_ds.front()};
                    recoil_ds.pop_front();
                    pair<int, int> key {confirmed_d->m_c, confirmed_d->m_d};
                    m_assigned_domains.push_back(key);
                }
                recoil_ds.push_back(cur_d);
                c_attempts_q.push_back(c_attempts);
                m_avail_cis.push_back(cur_avail_cis);
                if (m_regrow_ds.empty()) {
                    growth_q_empty = true;
                }
                else {
                    prev_p = cur_d->m_pos;
                    cur_d = m_regrow_ds.front();
                    m_regrow_ds.pop_front();
                    c_attempts = 0;
                    growthpoint = m_constraintpoints.is_growthpoint(cur_d);
                    if (growthpoint) {
                        cur_c_max_attempts = 1;
                        cur_avail_cis = {};
                        dir = m_regrow_ds.front()->m_d - cur_d->m_d;
                    }
                    else {
                        cur_c_max_attempts = m_max_c_attempts;
                        cur_avail_cis = m_all_cis;
                    }
                }
            }
            else {
                m_origami_system.unassign_domain(*cur_d);
                m_regrow_ds.push_front(cur_d);
                if (recoil_ds.empty()) {
                    recoil_q_empty = true;
                }
                else {
                    growthpoint = m_constraintpoints.is_growthpoint(cur_d);
                    cur_d = recoil_ds.front();
                    if (growthpoint) {
                        dir = m_regrow_ds.front()->m_d - cur_d->m_d;
                    }
                    recoil_ds.pop_front();
                    c_attempts = c_attempts_q.front();
                    c_attempts_q.pop_front();
                    cur_avail_cis = m_avail_cis.front();
                    m_avail_cis.pop_front();
                }
            }
        }
    }

    configT CTScaffoldRG::select_trial_config(
            Domain* d,
            vector<int>& avail_cis,
            bool growthpoint) {
        if (growthpoint) {
            int ci {m_random_gens.uniform_int(0, avail_cis.size() - 1)};
            configT c {m_all_configs[ci]};
            c.first = c.first + d->m_pos;
        }
        else {
        }

        return c;
    }

    bool CTScaffoldRG::test_config_open(configT config) {
        ;
    }

    void CTScaffoldRG::calc_rg_weights() {

        while (not m_regrow_ds.empty()) {
            int avail_cs {1};
            Domain* cur_d {m_regrow_ds.front()};
            m_regrow_ds.pop_front();
            int cur_c_attempts {c_attempts_q.front()};
            c_attempts_q.pop_front();
            while (cur_c_attempts != m_max_c_attempts) {
                configT config {select_trial_config(cur_d)};
                if (test_config_open(config)) {
                    bool c_avail {test_config_avail(config, m_regrow_ds)};
                    avail_cs += c_avail;
                }
            }
            m_weight *= avail_cs;
        }
    }

    bool CTScaffoldRG::test_config_avail(configT config) {
        bool c_avail;
        bool feeler_complete {false};
        deque<Domain*> feelers {};
        deque<int> c_attempts_q {};
        while (not feeler_complete) {
            Domain* cur_d {m_regrow_ds.front()};
            m_regrow_ds.pop_front();
            bool c_open {false};
            int c_attempts {0};
            while (not c_open and c_attempts != m_max_c_attempts) {
                c_attempts++;
                config = select_trial_config(cur_d);
                c_open = test_config_open(config); 
            }
            if (c_open) {
                int num_feelers {static_cast<int>(feelers.size())};
                if (num_feelers == m_max_recoils or m_regrow_ds.empty()) {
                    c_avail = true;
                    feeler_complete = true;
                }
                else {
                    cur_d = m_regrow_ds.front();
                    m_regrow_ds.pop_front();
                    c_attempts = 0;
                }
            }
            else {
                m_regrow_ds.push_back(cur_d);
                if (feelers.empty()) {
                    c_avail = false;
                    feeler_complete = true;
                }
                else {
                    cur_d = feelers.front();
                    c_attempts = c_attempts_q.front();
                    c_attempts_q.pop_front();
                }
            }
        }

        return c_avail;
    }
                    
    bool CTScaffoldRG::test_rg_acceptance() {
        long double ratio {m_new_weight / m_old_weight};
        bool accepted {false};
        if (test_acceptance(ratio)) {
            reset_origami();
            accepted = true;
        }

        // TODO change to match order of weight calculation
        else {
            m_modified_domains.clear();
            m_assigned_domains.clear();
            accepted = false;
        }

        return accepted;
    }
}
