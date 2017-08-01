// movetypes.h

#ifndef MOVETYPES_H
#define MOVETYPES_H

#include <iostream>
#include <iostream>
#include <unordered_map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "bias_functions.h"
#include "files.h"
#include "ideal_random_walk.h"
#include "domain.h"
#include "order_params.h"
#include "origami_system.h"
#include "parser.h"
#include "random_gens.h"
#include "utility.h"

namespace movetypes {

    using std::cout;
    using std::ostream;
    using std::unique_ptr;
    using std::pair;
    using std::set;
    using std::string;
    using std::unordered_map;
    using std::vector;

    using biasFunctions::SystemBiases;
    using domainContainer::Domain;
    using files::OrigamiOutputFile;
    using idealRandomWalk::IdealRandomWalks;
    using orderParams::SystemOrderParams;
    using origami::OrigamiSystem;
    using parser::InputParameters;
    using randomGen::RandomGens;
    using utility::VectorThree;

    struct MovetypeTracking {
        int attempts;
        int accepts;
    };

    class MCMovetype {
        public:
            MCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params);
            virtual ~MCMovetype() {};

            virtual bool attempt_move(long long int step) = 0;
            virtual void reset_origami();
            virtual void reset_internal();
            virtual void write_log_summary(ostream* log_stream) = 0;

            string get_label();
            int get_attempts();
            int get_accepts();

        protected:
            Domain* select_random_domain();
            int select_random_staple_identity();
            int select_random_staple_of_identity(int c_i_ident);
            VectorThree select_random_position(VectorThree p_prev);
            VectorThree select_random_orientation();
            bool test_acceptance(long double p_ratio);
            bool staple_is_connector(vector<Domain*> staple);
            set<int> find_staples(vector<Domain*> domains);
            bool scan_for_scaffold_domain(Domain*, set<int>& participating_chains);
            void write_config();

            virtual void add_external_bias() = 0;

            OrigamiSystem& m_origami_system;
            RandomGens& m_random_gens;
            IdealRandomWalks& m_ideal_random_walks;
            vector<OrigamiOutputFile*> m_config_files;
            string m_label;
            SystemOrderParams& m_ops;
            SystemBiases& m_biases;
            InputParameters& m_params;
            bool m_rejected {false};
            int m_config_output_freq;

            // Staple maxes
            int m_max_total_staples;
            int m_max_type_staples;

            // Lists of modified domains for move reversal
            vector<pair<int, int>> m_modified_domains {};
            vector<pair<int, int>> m_assigned_domains {};
            vector<int> m_added_chains {};

            // Domains mapped to previous configs
            unordered_map<pair<int, int>, VectorThree> m_prev_pos {};
            unordered_map<pair<int, int>, VectorThree> m_prev_ore {};

            // Modifier to correct for excluded volume and overcounting
            double m_modifier {1};

            MovetypeTracking m_general_tracker {0, 0};
            long long int m_step {0};
    };

    class IdentityMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move(long long int) {return true;};
    };

    class RegrowthMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;

        protected:
            virtual double set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            void grow_staple(int d_i_index, vector<Domain*> selected_chain);
            virtual void grow_chain(vector<Domain*> domains) = 0;
            pair<Domain*, Domain*> select_new_growthpoint(vector<Domain*> selected_chain);
    };

    template <typename T>
    void add_tracker(
            T tracker,
            unordered_map<T, MovetypeTracking>& trackers,
            bool accepted) {
        if (trackers.find(tracker) == trackers.end()) {
            trackers[tracker] = {1, accepted};
        }
        else {
            trackers[tracker].attempts++;
            trackers[tracker].accepts += accepted;
        }
    }
}

#endif // MOVETYPES_H
