// movetypes.h

#ifndef MOVETYPES_H
#define MOVETYPES_H

#include <memory>

#include "origami_system.h"

using namespace Origami;
using std::unique_ptr;

namespace Movetypes {

    class MCMovetype {
        public:
            MCMovetype(OrigamiSystem& origami_system) : m_origami_system
                    {origami_system} {};

            virtual bool attempt_move() = 0;
            void reset_origami();

        protected:
            OrigamiSystem& m_origami_system;

            // Lists of modified domains for move reversal
            vector<Domain*> m_modified_domains {};
            vector<Domain*> m_assigned_domains {};
            vector<int> m_added_chains {};
            vector<Domain*> m_added_domains {};

            // Contains the chain index and the identity, respectively
            vector<pair<int, int>> m_deleted_chains {};
            vector<Domain*> m_deleted_domains {};

            // Domains mapped to previous configs
            unordered_map<Domain*, VectorThree> m_prev_pos {};
            unordered_map<Domain*, VectorThree> m_prev_ore {};

            // Modifier to correct for excluded volume and overcounting
            double m_modifier {1};

            // Shared methods
            Domain* select_random_domain();
            int select_random_staple_identity();
            int select_random_staple_of_identity(int c_i_ident);
            double set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            void calc_overcount(Domain& domain);
    };

    template<typename T>
    unique_ptr<MCMovetype> movetype_constructor(OrigamiSystem origami_system);

    using MovetypeConstructor = unique_ptr<MCMovetype> (*)(OrigamiSystem origami_system);

    struct Movetype {
        MovetypeConstructor identity {movetype_constructor<MCMovetype>};
        MovetypeConstructor orientation_rotation {movetype_constructor<MCMovetype>};
        MovetypeConstructor staple_exchange {movetype_constructor<MCMovetype>};
        MovetypeConstructor cb_staple_regrowth {movetype_constructor<MCMovetype>};
        MovetypeConstructor ctcb_scaffold_regrowth {movetype_constructor<MCMovetype>};
    };

    bool test_acceptance(double p_ratio, double modifier);

    VectorThree select_random_position(VectorThree p_prev);
    VectorThree select_random_orientation();

    class IdentityMCMovetype: public MCMovetype {
        public:
            bool attempt_move() {return true;};
    };

    class OrientationRotationMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move();
    };

    class MetMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
        protected:
            double m_delta_e {0};

            void grow_staple(int d_i_index, vector<Domain*> selected_chain);
            void grow_chain(vector<Domain*> domains);
    };

    class StapleExchangeMCMovetype: public MetMCMovetype {
        public:
            using MetMCMovetype::MetMCMovetype;
            bool attempt_move();

        protected:

            // These can be overidden for a derived class the excludes misbinding
            int preconstrained_df {0};
            int m_insertion_sites {m_origami_system.m_num_domains};

        private:
            bool staple_insertion_accepted(int c_i_ident);
            bool staple_deletion_accepted(int c_i_ident);
            bool insert_staple();
            bool delete_staple();
            void select_and_set_growth_point(Domain* growth_domain_new);
    };

    class CBStapleRegrowthMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move();
    };

    class CTCBScaffoldRegrowth: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
    };

}

#endif // MOVETYPES_H
