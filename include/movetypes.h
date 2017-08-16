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
#include "top_constraint_points.h"
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
    using topConstraintPoints::Constraintpoints;
    using utility::VectorThree;

    /**
      * Tracker for general movetype info
      */
    struct MovetypeTracking {
        int attempts;
        int accepts;
    };

    /**
      * Parent class of all movetypes
      * 
      * It serves both to define the interface and provides some shared
      * implementation details, including shared data.
      */
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

            /** Attempt move and return result */
            bool attempt_move(long long int step);

            /** Reset origami system to state before move attempt */
            virtual void reset_origami();

            /** Write summary of movetype attempts to file */
            virtual void write_log_summary(ostream* log_stream) = 0;

            string get_label();
            int get_attempts();
            int get_accepts();

        protected:
            // Methods requiring definition for all movetypes
            virtual void add_external_bias() = 0;
            virtual bool internal_attempt_move() = 0;
            virtual void add_tracker(bool accepted) = 0;

            // Methods that probably need to be overriden
            virtual void reset_internal();

            // Shared general utility functions

            /** Return a random domain in the system
              *
              * Uniform over the set of all domains, rather than uniform over
              * chains and then over the constituent domains
              */
            Domain* select_random_domain();

            /** Return a random staple type */
            int select_random_staple_identity();

            /** Return index to a random staple staple of given type */
            int select_random_staple_of_identity(int c_i_ident);

            /** Return a random neighbour lattice site position to given */
            VectorThree select_random_position(VectorThree p_prev);

            /** Return a random unit orientation vector */
            VectorThree select_random_orientation();

            /** Test if move accepted with given probality */
            bool test_acceptance(long double p_ratio);

            /** Test if staple is anchoring other staples to the system
              *
              * Connecting staples, or anchor staples, are staples that if
              * removed would leave the system in a state with staples
              * unconnected to the scaffold
              */
            bool staple_is_connector(vector<Domain*> staple);

            set<int> find_staples(vector<Domain*> domains);

            /** Write the config to move file
              *
              * The move file is for recording the configuration as the chains
              * are grown.
              */
            void write_config();

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

        private:
            bool scan_for_scaffold_domain(Domain*, set<int>& participating_chains);
    };

    /** For debugging purposes */
    class IdentityMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move(long long int) {return true;};
    };

    /**
      * Parent class of moves with chain regrowth
      */
    class RegrowthMCMovetype:
        virtual public MCMovetype {

        public:
            using MCMovetype::MCMovetype;

        protected:
            virtual void grow_chain(vector<Domain*> domains) = 0;

            double set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            void grow_staple(int d_i_index, vector<Domain*> selected_chain);
            pair<Domain*, Domain*> select_new_growthpoint(vector<Domain*> selected_chain);
    };

    /**
      * Parent class of moves with constant topology chain regrowth
      */
    class CTRegrowthMCMovetype:
        virtual public RegrowthMCMovetype {

        public:
            CTRegrowthMCMovetype(int num_excluded_staples);

        protected:
            void sel_excluded_staples();

            /** Select scaffold segment to be regrown */
            vector<Domain*> select_indices(vector<Domain*>);

            /** Return segment of linear domains given endpoints */
            vector<Domain*> select_cyclic_segment(
                    Domain* start_domain,
                    Domain* end_domain);

            /**
              * Return segment of cyclic domains given endpoints
              *
              * Selects a direction of regrowth with uniform probability
              */
            vector<Domain*> select_linear_segment(
                    Domain* start_domain,
                    Domain* end_domain);

            /** Select endpoints of a chain segment of given length */
            pair<int, int> select_endpoints(
                    const int array_size,
                    const int mar, // Number of elements to exclude at ends
                    const int min_size); // Minimum size of selection (excluding endpoints)

            /** Check if excluded staples are still bound to system */
            bool excluded_staples_bound();

            int m_dir;
            int m_num_excluded_staples;
            vector<int> m_excluded_staples;
            Constraintpoints m_constraintpoints {m_origami_system,
                    m_ideal_random_walks}; // For fixed end biases
            vector<Domain*> m_scaffold {};
    };

    /**
      * Template function for updating movetype specific trackers
      */
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
