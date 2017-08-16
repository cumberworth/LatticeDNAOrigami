// rg_movetypes.h

#ifndef RG_MOVETYPES_H
#define RG_MOVETYPES_H

#include <deque>

#include "movetypes.h"

namespace movetypes {

    using std::deque;

    typedef pair<VectorThree, VectorThree> configT;

    /**
     * Recoil regrowth of scaffold segment and bound staples
     */
    class CTScaffoldRG:
        public CTRegrowthMCMovetype {

        public:
            CTScaffoldRG(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params,
                    int num_excluded_staples,
                    int max_num_recoils);

            void write_log_summary(ostream* log_entry) override;
            void reset_internal() override;

        private:
            bool internal_attempt_move() override;
            void add_tracker(bool accepted) override;

            void unassign_domains();

            /** Update the external bias without adding */
            void update_external_bias();
            void recoil_regrow();
            configT select_trial_config(
                    Domain* d,
                    vector<int>& avail_cis,
                    bool growthpoint);
            bool test_config_open(configT config);
            void calc_rg_weights();
            bool test_config_avail(configT config);
            bool test_rg_acceptance();

            const vector<configT> m_all_configs;
            const vector<int> m_all_cis;

            // Movetype parameters
            int m_max_recoils;
            int m_max_c_attempts;

            // Move states (reset for each move)
            vector<Domain*> m_sel_scaf_doms {};
            deque<Domain*> m_regrow_ds {};
            deque<int> m_c_attempts_q {};
            deque<vector<int>> m_avail_cis {};
    };
    
}

#endif // RG_MOVETYPES_H
