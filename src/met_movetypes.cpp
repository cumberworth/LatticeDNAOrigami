// met_movetypes.cpp

#include <map>
#include <set>
#include <utility>

#include "met_movetypes.h"

namespace MetMovetypes {

    using std::map;
    using std::pair;
    using std::set;

    using Movetypes::add_tracker;
    using Utility::VectorThree;

    void MetMCMovetype::reset_internal() {
        MCMovetype::reset_internal();
        m_delta_e = 0;
    }

    void MetMCMovetype::grow_chain(vector<Domain*> domains) {
        for (size_t i {1}; i != domains.size(); i++) {
            Domain* domain {domains[i]};
            Domain* prev_domain {domains[i - 1]};
            VectorThree new_p {select_random_position(prev_domain->m_pos)};
            VectorThree new_o {select_random_orientation()};
            m_delta_e += m_origami_system.set_domain_config(*domain, new_p, new_o);
            if (m_origami_system.m_constraints_violated) {
                m_rejected = true;
                break;
            }
            else {
                pair<int, int> key {domain->m_c, domain->m_d};
                m_assigned_domains.push_back(key);
                write_config();
            }
        }
    }

    void MetMCMovetype::unassign_domains(vector<Domain*> domains) {
        for (auto domain: domains) {
            pair<int, int> key {domain->m_c, domain->m_d};
            m_prev_pos[key] = domain->m_pos;
            m_prev_ore[key] = domain->m_ore;
            m_modified_domains.push_back(key);
            m_delta_e += m_origami_system.unassign_domain(*domain);
        }
    }

    void MetMCMovetype::add_external_bias() {
        double total_bias {m_system_bias.calc_bias()};
        m_delta_e += total_bias;

        return;
    }

    void MetMCMovetype::subtract_external_bias() {
        double total_bias {m_system_bias.calc_bias()};
        m_delta_e -= total_bias;

        return;
    }

    bool MetStapleExchangeMCMovetype::attempt_move(long long int step) {
        m_step = step;
        write_config();
        m_general_tracker.attempts++;
        subtract_external_bias();

        // Select inertion or deletion with equal frequency
        bool accept;
        if (m_random_gens.uniform_real() < 0.5) {
            m_tracker.staple_insertion = true;
            accept = insert_staple();
        }
        else {
            m_tracker.staple_insertion = false;
            accept = delete_staple();
        }
        
        return accept;
    }

    void MetStapleExchangeMCMovetype::write_log_summary(ostream* log_stream) {

        // Insertion of each staple type
        map<int, int> insertion_attempts {};
        map<int, int> insertion_accepts {};
        map<int, int> deletion_attempts {};
        map<int, int> deletion_accepts {};
        set<int> staple_types {};
        for (auto tracker: m_tracking) {
            auto info = tracker.first;
            auto counts = tracker.second;
            staple_types.insert(info.staple_type);
            if (info.staple_insertion) {
                insertion_attempts[info.staple_type] = counts.attempts;
                insertion_accepts[info.staple_type] = counts.accepts;
            }
            else if (not info.no_staples) {
                deletion_attempts[info.staple_type] = counts.attempts;
                deletion_accepts[info.staple_type] = counts.accepts;
            }
        }

        for (auto st: staple_types) {
            *log_stream << "    Staple type: " << st << "\n";
            int iats {insertion_attempts[st]};
            int iacs {insertion_accepts[st]};
            float ifreq {static_cast<float>(iacs) / iats};
            *log_stream << "        Insertion attempts: " << iats << "\n";
            *log_stream << "        Insertion accepts: " << iacs << "\n";
            *log_stream << "        Insertion frequency: " << ifreq << "\n";
            int dats {deletion_attempts[st]};
            int dacs {deletion_accepts[st]};
            float dfreq {static_cast<float>(dacs) / dats};
            *log_stream << "        Deletion attempts: " << dats << "\n";
            *log_stream << "        Deletion accepts: " << dacs << "\n";
            *log_stream << "        Deletion frequency: " << dfreq << "\n";
        }
    }

    void MetStapleExchangeMCMovetype::reset_internal() {
        MetMCMovetype::reset_internal();
        m_insertion_sites = m_origami_system.num_domains();
    }

    bool MetStapleExchangeMCMovetype::staple_insertion_accepted(int c_i_ident) {
        add_external_bias();
        double boltz_factor {exp(-m_delta_e)};
        int Ni_new {m_origami_system.num_staples_of_ident(c_i_ident)};

        // Correct for extra states from additional staple domains
        size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
        int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
        double extra_states {pow(6, extra_df)};
        double ratio {extra_states / Ni_new * boltz_factor};

        // Correct for insertion into subset of volume
        m_modifier *= m_insertion_sites / m_origami_system.m_volume;

        // Correct for considering only 1 of staple length ways insertion could occur
        m_modifier *= staple_length;

        // Exchange probability multiplier
        m_modifier *= m_params.m_exchange_mult;

        return test_acceptance(ratio);
    }

    bool MetStapleExchangeMCMovetype::staple_deletion_accepted(int c_i_ident) {
        // SUPER HACK
        double total_bias {m_system_bias.calc_bias_delete_case()};
        m_delta_e += total_bias;

        double boltz_factor {exp(-m_delta_e)};
        int Ni {m_origami_system.num_staples_of_ident(c_i_ident)};

        // Correct for extra states from additional staple domains
        size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
        double extra_df {2 * static_cast<double>(staple_length) - 1 - preconstrained_df};
        double extra_states {pow(6, extra_df)};
        double ratio {Ni / extra_states * boltz_factor};

        // Exchange probability multiplier
        m_modifier *= m_params.m_exchange_mult;

        return test_acceptance(ratio);
    }

    bool MetStapleExchangeMCMovetype::insert_staple() {
        bool accepted {false};

        // Select and add chain of random identity
        int c_i_ident {select_random_staple_identity()};
        m_tracker.staple_type = c_i_ident;

        // Check if number staples exceeds max allowed
        if (m_origami_system.num_staples() == m_max_total_staples) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }
        if (m_origami_system.num_staples_of_ident(c_i_ident) == m_max_type_staples) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        int c_i {m_origami_system.add_chain(c_i_ident)};
        m_added_chains.push_back(c_i);

        // Assume that add_chain always adds to end of m_domains
        vector<Domain*> selected_chain {m_origami_system.get_last_chain()};

        // Select growth points on chains
        pair<Domain*, Domain*> growthpoint {select_new_growthpoint(selected_chain)};

        m_delta_e += set_growth_point(*growthpoint.first, *growthpoint.second);
        if (m_rejected) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        grow_staple(growthpoint.first->m_d, selected_chain);
        if (m_rejected) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        accepted = staple_insertion_accepted(c_i_ident);
        add_tracker(m_tracker, m_tracking, accepted);
        m_general_tracker.accepts += accepted;
        return accepted;
    }

    bool MetStapleExchangeMCMovetype::delete_staple() {
        bool accepted {false};

        // Select random identity and staple of that type
        int c_i_ident {select_random_staple_identity()};
        m_tracker.staple_type = c_i_ident;
        int c_i {select_random_staple_of_identity(c_i_ident)};
        if (m_rejected) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        // Reject if staple is connector
        vector<Domain*> staple {m_origami_system.get_chain(c_i)};
        if (staple_is_connector(staple)) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        // Unassign domains and test acceptance
        unassign_domains(staple);
        m_delta_e += m_origami_system.check_delete_chain(c_i);
        accepted = staple_deletion_accepted(c_i_ident);
        add_tracker(m_tracker, m_tracking, accepted);
        m_general_tracker.accepts += accepted;
        if (accepted) {
            m_origami_system.delete_chain(c_i);
        }

        return accepted;
    }

    bool MetStapleRegrowthMCMovetype::attempt_move(long long int step) {
        m_step = step;
        write_config();
        bool accepted {false};
        m_general_tracker.attempts++;
        subtract_external_bias();

        // No staples to regrow
        if (m_origami_system.num_staples() == 0) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
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

        unassign_domains(selected_chain);

        // Grow staple
        m_delta_e += set_growth_point(*growthpoint.first, *growthpoint.second);
        if (m_rejected) {
            accepted = false;
            return accepted;
        }
        grow_staple(growthpoint.first->m_d, selected_chain);
        if (m_rejected) {
            add_tracker(m_tracker, m_tracking, accepted);
            return accepted;
        }

        add_external_bias();
        double boltz_factor {exp(-(m_delta_e))};
        accepted = test_acceptance(boltz_factor);
        m_general_tracker.accepts++;
        add_tracker(m_tracker, m_tracking, accepted);
        return accepted;
    }

    void MetStapleRegrowthMCMovetype::write_log_summary(ostream* log_stream) {

        // Insertion of each staple type
        map<int, int> attempts {};
        map<int, int> accepts {};
        set<int> staple_types {};
        for (auto tracker: m_tracking) {
            auto info = tracker.first;
            auto counts = tracker.second;
            staple_types.insert(info.staple_type);
            if (not info.no_staples) {
                attempts[info.staple_type] = counts.attempts;
                accepts[info.staple_type] = counts.accepts;
            }
        }

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
}
