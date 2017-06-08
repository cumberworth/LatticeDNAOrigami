// movetypes.h

#ifndef MOVETYPES_H
#define MOVETYPES_H

#include <iostream>
#include <set>

#include <memory>

#include "parser.h"
#include "random_gens.h"
#include "origami_system.h"
#include "ideal_random_walk.h"
#include "order_params.h"

using std::cout;
using std::unique_ptr;
using std::set;

using namespace Parser;
using namespace Origami;
using namespace RandomGen;
using namespace IdealRandomWalk;
using namespace OrderParams;

namespace Movetypes {

    // Movetype classes

    class MCMovetype {
        public:
            MCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    InputParameters& params) :
                    m_origami_system {origami_system},
                    m_random_gens {random_gens},
                    m_ideal_random_walks {ideal_random_walks},
                    m_params {params},
                    //HACK
                    m_system_bias {*origami_system.get_system_biases()},
                    m_max_total_staples {params.m_max_total_staples},
                    m_max_type_staples {params.m_max_type_staples} {}
            virtual ~MCMovetype() {};

            virtual bool attempt_move() = 0;
            virtual void reset_origami();
            virtual void reset_internal();

            virtual string m_label() {return "MCMovetype";};

        protected:
            OrigamiSystem& m_origami_system;
            RandomGens& m_random_gens;
            IdealRandomWalks& m_ideal_random_walks;
            InputParameters& m_params;
            bool m_rejected {false};
            //HACK
            SystemBiases& m_system_bias;

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

            // Shared methods
            Domain* select_random_domain();
            int select_random_staple_identity();
            int select_random_staple_of_identity(int c_i_ident);
            VectorThree select_random_position(VectorThree p_prev);
            VectorThree select_random_orientation();
            bool test_acceptance(long double p_ratio);
            bool staple_is_connector(vector<Domain*> staple);
            bool scan_for_scaffold_domain(Domain*, set<int>& participating_chains);

            virtual void update_bias(int sign) = 0;
    };

    class IdentityMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move() {return true;};

            string m_label() {return "IdentityMCMovetype";};
    };

    class RegrowthMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;

            string m_label() {return "RegrowthMCMovetype";};
        protected:
            virtual double set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            void grow_staple(int d_i_index, vector<Domain*> selected_chain);
            virtual void grow_chain(vector<Domain*> domains) = 0;
            pair<Domain*, Domain*> select_new_growthpoint(vector<Domain*> selected_chain);
    };
}

#endif // MOVETYPES_H
