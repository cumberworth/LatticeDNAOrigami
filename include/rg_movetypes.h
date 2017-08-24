// rg_movetypes.h

#ifndef RG_MOVETYPES_H
#define RG_MOVETYPES_H

#include <list>

#include "movetypes.h"

namespace movetypes {

    using std::deque;
    using std::list;

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
                    int max_num_recoils,
                    int max_c_attempts);

            void write_log_summary(ostream* log_entry) override;
            void reset_internal() override;

        private:
            void add_external_bias() override;
            bool internal_attempt_move() override;
            void add_tracker(bool accepted) override;

            void unassign_and_save_domains();
            void unassign_domains();

            /** Update the external bias without adding */
            void update_external_bias();
            void setup_constraints();
            void setup_for_calc_new_weights();
            void setup_for_calc_old_weights();
            void recoil_regrow();
            void set_config(Domain* d, configT c);
            void prepare_for_growth();
            void prepare_for_regrowth();
            void restore_endpoints();
            configT select_trial_config();
            double calc_p_config_open(configT c);
            bool test_config_open(double p_c_open);
            void calc_weights();
            bool test_config_avail();
            void calc_old_c_opens();
            bool test_rg_acceptance();

            const vector<configT> m_all_configs;
            const unordered_map<configT, int> m_config_to_i;
            list<int> m_all_cis;

            // Movetype parameters
            int m_max_recoils;
            int m_max_c_attempts;

            // Overall move states (reset for each move)
            vector<Domain*> m_sel_scaf_doms {}; // Selected scaffold region
            vector<Domain*> m_regrow_ds {}; // Domains to regrow
            vector<int> m_c_attempts_q {}; // Num. attempted configs per domain
            vector<int> m_c_attempts_wq {}; // Store of above for weight calc
            vector<list<int>> m_avail_cis_q {}; // Available configs per domain
            vector<list<int>> m_avail_cis_wq {}; // Store of above for weight calc
            vector<double> m_c_opens; // Configuration open probability
            unordered_map<pair<int, int>, VectorThree> m_old_pos {};
            unordered_map<pair<int, int>, VectorThree> m_old_ore {};
            vector<vector<VectorThree>> m_erased_endpoints_q {};

            double m_delta_e; // Energy change
            double m_delta_e_new; // Energy change of new config
            double m_weight; // RG weight
            double m_weight_new; // RG weight of new config

            // States for domain currently being regrown
            unsigned int m_di; // Index into m_regrow_ds of the current domain
            Domain* m_d;  // Domain to be regrown
            Domain* m_ref_d; // Reference domain
            bool m_prev_gp; // If the previous domain is a growthpoint
            int m_c_attempts; // Number of configs tried for current domain
            int m_d_max_c_attempts; // Max number of configs to be tried for current domain
            list<int> m_avail_cis; // Available configurations

    };
    
}

#endif // RG_MOVETYPES_H
