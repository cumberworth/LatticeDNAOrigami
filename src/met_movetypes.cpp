// met_movetypes.cpp

#include <map>
#include <set>
#include <utility>

#include "met_movetypes.h"

namespace movetypes {

    using std::map;
    using std::pair;
    using std::set;
    using utility::SimulationMisuse;
    using utility::Occupancy;

	MetMCMovetype::MetMCMovetype(
                OrigamiSystem& origami_system,
                RandomGens& random_gens,
                IdealRandomWalks& ideal_random_walks,
                vector<OrigamiOutputFile*> config_files,
                string label,
                SystemOrderParams& ops,
                SystemBiases& biases,
                InputParameters& params) :
        MCMovetype(origami_system, random_gens, ideal_random_walks,
                config_files, label, ops, biases, params),
        RegrowthMCMovetype(origami_system, random_gens, ideal_random_walks,
                config_files, label, ops, biases, params) {
	}

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

    void MetMCMovetype::add_external_bias() {
        m_ops.update_move_params();
        double ex_bias {m_biases.calc_move()};
        m_delta_e += ex_bias;

        return;
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

    MetStapleExchangeMCMovetype::MetStapleExchangeMCMovetype(
            OrigamiSystem& origami_system,
            RandomGens& random_gens,
            IdealRandomWalks& ideal_random_walks,
            vector<OrigamiOutputFile*> config_files,
            string label,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params,
            vector<double> exchange_mults,
            bool adaptive_exchange):
            MCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            RegrowthMCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            MetMCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            m_adaptive_exchange {adaptive_exchange},
            m_allow_nonsensical_ps {m_params.m_allow_nonsensical_ps},
            m_exchange_mults {exchange_mults} {

        size_t num_missing_exchange_mults {
                m_origami_system.m_identities.size() -
                m_exchange_mults.size() - 1};
        if (num_missing_exchange_mults > 0) {
            for (size_t i {0}; i != num_missing_exchange_mults; i++) {
                m_exchange_mults.push_back(1);
            }
        }
    }

    void MetStapleExchangeMCMovetype::reset_internal() {
        MetMCMovetype::reset_internal();
        m_insertion_sites = m_origami_system.num_domains();
        m_staple_bound = false;
    }

    void MetStapleExchangeMCMovetype::write_log_summary(ostream* log_stream) {
        write_log_summary_header(log_stream);

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

        *log_stream << "       Exchange multipliers\n        ";
        for (size_t i {0}; i != staple_types.size(); i++) {
            *log_stream << m_exchange_mults[i];
            *log_stream << ", ";
        }
        *log_stream << "\n";

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

    bool MetStapleExchangeMCMovetype::internal_attempt_move() {
        bool accepted {false};

        // Select inertion or deletion with equal frequency
        if (m_random_gens.uniform_real() < 0.5) {
            m_tracker.staple_insertion = true;
            accepted = insert_staple();
        }
        else {
            m_tracker.staple_insertion = false;
            accepted = delete_staple();
        }
        
        return accepted;
    }

    void MetStapleExchangeMCMovetype::add_tracker(bool accepted) {
        movetypes::add_tracker(m_tracker, m_tracking, accepted);
    }

    bool MetStapleExchangeMCMovetype::staple_insertion_accepted(int c_i_ident,
            int num_staple_bd) {

        add_external_bias();
        double boltz_factor {exp(-m_delta_e)};
        int Ni_new {m_origami_system.num_staples_of_ident(c_i_ident)};

        // Correct for extra states from additional staple domains
        size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
        //int extra_df {2 * static_cast<int>(staple_length) - 1 - preconstrained_df};
        //double extra_states {staple_length * pow(6, extra_df)};
        double pratio {staple_length / Ni_new * boltz_factor};

        // Correct for insertion into subset of volume
        pratio *= m_insertion_sites / m_origami_system.m_volume;

        // Correct for overcounting multiply bound staples
        pratio /= num_staple_bd;

        // Exchange probability multiplier
        bool accepted;
        if (m_staple_bound) {
            m_modifier *= m_exchange_mults[c_i_ident - 1];

            // Check if nonsensical probabilities will result
            if (m_modifier * fmin(1, pratio) > 1) {
                if (m_adaptive_exchange) {
                    m_exchange_mults[c_i_ident - 1] /= 10;
                    accepted = false;
                }
                else if (m_allow_nonsensical_ps) {
                    accepted = true;
                }
                else {
                    cout << "Nonsensical exchange probability detected\n";
                    throw SimulationMisuse {};
                }
            }
            else {
                accepted = test_acceptance(pratio);
            }
        }
        else {
            accepted = test_acceptance(pratio);
        }

        return accepted;
    }

    bool MetStapleExchangeMCMovetype::staple_deletion_accepted(int c_i_ident,
            int num_staple_bd) {

        // This is ugly but necessary as I don't remove staple unles accepted
        m_origami_system.temp_reduce_staples_by_one();
        add_external_bias();
        m_origami_system.undo_reduce_staples_by_one();
        double boltz_factor {exp(-m_delta_e)};
        int Ni {m_origami_system.num_staples_of_ident(c_i_ident)};

        // Correct for extra states from additional staple domains
        size_t staple_length {m_origami_system.m_identities[c_i_ident].size()};
        //double extra_df {2 * static_cast<double>(staple_length) - 1 - preconstrained_df};
        //double extra_states {staple_length * pow(6, extra_df)};
        double pratio {Ni / staple_length * boltz_factor};

        // Correct for insertion into subset of volume
        pratio *= m_origami_system.m_volume / (m_insertion_sites -
                staple_length);

        // Correct for overcounting multiply bound staples
        pratio *= num_staple_bd;

        // Exchange probability multiplier
        bool accepted;
        if (m_staple_bound) {
            m_modifier *= m_exchange_mults[c_i_ident - 1];

            // Check if nonsensical probabilities will result
            if (m_modifier * fmin(1, pratio) > 1) {
                if (m_adaptive_exchange) {
                    cout << c_i_ident << m_modifier << " " << pratio << "\n";
                    m_exchange_mults[c_i_ident - 1] /= 10;
                    accepted = false;
                }
                else if (m_allow_nonsensical_ps) {
                    accepted = true;
                }
                else {
                    cout << "Nonsensical exchange probability detected\n";
                    throw SimulationMisuse {};
                }
            }
            else {
                accepted = test_acceptance(pratio);
            }
        }
        else {
            accepted = test_acceptance(pratio);
        }

        return accepted;
    }

    bool MetStapleExchangeMCMovetype::insert_staple() {
        bool accepted {false};

        // Select and add chain of random identity
        int c_i_ident {select_random_staple_identity()};
        m_tracker.staple_type = c_i_ident;

        // Check if number staples exceeds max allowed
        if (m_origami_system.num_staples() == m_max_total_staples) {
            return accepted;
        }
        if (m_origami_system.num_staples_of_ident(c_i_ident) == m_max_type_staples) {
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
            return accepted;
        }

        grow_staple(growthpoint.first->m_d, selected_chain);
        if (m_rejected) {
            return accepted;
        }

        vector<Domain*> staple {m_origami_system.get_chain(c_i)};
        for (size_t i {0}; i != staple.size(); i++) {
            Domain* d {staple[i]};
            if (d->m_state == Occupancy::bound) {
                m_staple_bound = true;
                break;
            }
        }

        accepted = staple_insertion_accepted(c_i_ident,
                num_bound_staple_domains(staple));
        return accepted;
    }

    bool MetStapleExchangeMCMovetype::delete_staple() {
        bool accepted {false};

        // Select random identity and staple of that type
        int c_i_ident {select_random_staple_identity()};
        m_tracker.staple_type = c_i_ident;
        int c_i {select_random_staple_of_identity(c_i_ident)};
        if (m_rejected) {
            return accepted;
        }

        // Reject if staple is connector
        vector<Domain*> staple {m_origami_system.get_chain(c_i)};
        if (staple_is_connector(staple)) {
            return accepted;
        }
        for (size_t i {0}; i != staple.size(); i++) {
            Domain* d {staple[i]};
            if (d->m_state == Occupancy::bound) {
                m_staple_bound = true;
                break;
            }
        }

        // Unassign domains and test acceptance
        int num_staple_bd {num_bound_staple_domains(staple)};
        unassign_domains(staple);
        accepted = staple_deletion_accepted(c_i_ident,
                num_staple_bd);
        if (accepted) {
            m_origami_system.delete_chain(c_i);
        }

        return accepted;
    }

    MetStapleRegrowthMCMovetype::MetStapleRegrowthMCMovetype(
            OrigamiSystem& origami_system,
            RandomGens& random_gens,
            IdealRandomWalks& ideal_random_walks,
            vector<OrigamiOutputFile*> config_files,
            string label,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params):
            MCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            RegrowthMCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params),
            MetMCMovetype(origami_system, random_gens, ideal_random_walks,
                    config_files, label, ops, biases, params) {
    }

    void MetStapleRegrowthMCMovetype::write_log_summary(ostream* log_stream) {
        write_log_summary_header(log_stream);

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

    bool MetStapleRegrowthMCMovetype::internal_attempt_move() {
        bool accepted {false};

        // No staples to regrow
        if (m_origami_system.num_staples() == 0) {
            m_tracker.no_staples = true;
            return accepted;
        }

        // Select a staple to regrow
        int c_i_index {m_random_gens.uniform_int(1,
                m_origami_system.num_staples())};
        vector<Domain*> selected_chain {
                m_origami_system.get_chains()[c_i_index]};
        m_tracker.staple_type = selected_chain[0]->m_c_ident;

        // Reject if staple is connector
        if (staple_is_connector(selected_chain)) {
            return accepted;
        }

        // Select growth points on chains
        auto bound_domains = find_bound_domains(selected_chain);
        domainPairT growthpoint {select_old_growthpoint(bound_domains)};

        unassign_domains(selected_chain);

        // Grow staple
        m_delta_e += set_growth_point(*growthpoint.first, *growthpoint.second);
        if (m_rejected) {
            accepted = false;
            return accepted;
        }
        grow_staple(growthpoint.first->m_d, selected_chain);
        if (m_rejected) {
            return accepted;
        }

        add_external_bias();
        int new_num_staple_bd {num_bound_staple_domains(selected_chain)};
        double boltz_factor {exp(-(m_delta_e))};
        double pratio {boltz_factor * bound_domains.size() / new_num_staple_bd};
        accepted = test_acceptance(pratio);
        return accepted;
    }

    void MetStapleRegrowthMCMovetype::add_tracker(bool accepted) {
        movetypes::add_tracker(m_tracker, m_tracking, accepted);
    }
}
