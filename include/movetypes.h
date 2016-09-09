// movetypes.h

#ifndef MOVETYPES_H
#define MOVETYPES_H

#include <iostream>

#include <memory>

#include "random_gens.h"
#include "origami_system.h"

using std::cout;
using std::unique_ptr;

using namespace Origami;
using namespace RandomGen;

namespace Movetypes {

    // Movetype classes

    class MCMovetype {
        public:
            MCMovetype(OrigamiSystem& origami_system, RandomGens& random_gens) : m_origami_system
                    {origami_system}, m_random_gens {random_gens} {};

            virtual bool attempt_move() = 0;
            void reset_origami();

        protected:
            OrigamiSystem& m_origami_system;
            RandomGens& m_random_gens;

            // Lists of modified domains for move reversal
            vector<pair<int, int>> m_modified_domains {};
            vector<pair<int, int>> m_assigned_domains {};
            vector<int> m_added_chains {};

            // Contains the chain index and the identity, respectively
            vector<pair<int, int>> m_deleted_chains {};
            map<int, vector<int>> m_deleted_domains {};

            // Domains mapped to previous configs
            unordered_map<pair<int, int>, VectorThree> m_prev_pos {};
            unordered_map<pair<int, int>, VectorThree> m_prev_ore {};

            // Modifier to correct for excluded volume and overcounting
            double m_modifier {1};

            // Shared methods
            VectorThree select_random_position(VectorThree p_prev);
            VectorThree select_random_orientation();
            bool test_acceptance(double p_ratio);
            Domain* select_random_domain();
            int select_random_staple_identity();
            int select_random_staple_of_identity(int c_i_ident);
    };

    class IdentityMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move() {return true;};
    };

    class OrientationRotationMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move();
    };

    class RegrowthMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
        protected:
            double set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            void grow_staple(int d_i_index, vector<Domain*> selected_chain);
            virtual void grow_chain(vector<Domain*> domains) = 0;
    };

    class MetMCMovetype: public RegrowthMCMovetype {
        public:
            using RegrowthMCMovetype::RegrowthMCMovetype;
        protected:
            double m_delta_e {0};

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

    class CBMCMovetype: public RegrowthMCMovetype {
        public:
            using RegrowthMCMovetype::RegrowthMCMovetype;
        protected:
            double m_bias {1};
            bool m_regrow_old {false};
            unordered_map<pair<int, int>, VectorThree> m_old_pos {};
            unordered_map<pair<int, int>, VectorThree> m_old_ore {};

    };

    class CBStapleRegrowthMCMovetype: public CBMCMovetype {
        public:
            using CBMCMovetype::CBMCMovetype;
            bool attempt_move();
        private:
            void grow_chain(vector<Domain*> domains);
            void calculate_biases(
                    Domain& domain,
                    VectorThree p_prev,
                    vector<pair<VectorThree, VectorThree>>& configs,
                    vector<double>& bfactors);
            vector<double> calc_bias(vector<double> bfactors);
            void select_and_set_new_config(
                    Domain& domain,
                    vector<double> weights,
                    vector<pair<VectorThree, VectorThree>> configs);
            void select_and_set_old_config(Domain& domain);
    };

    class CTCBScaffoldRegrowth: public CBMCMovetype {
        public:
            using CBMCMovetype::CBMCMovetype;
    };

    // Movetype construction

    template<typename T>
    unique_ptr<MCMovetype> movetype_constructor(OrigamiSystem& origami_system,
            RandomGens& random_gens);

    using MovetypeConstructor = unique_ptr<MCMovetype> (*)(OrigamiSystem&
            origami_system, RandomGens& random_gens);

    struct Movetype {
        MovetypeConstructor identity {movetype_constructor<IdentityMCMovetype>};
        MovetypeConstructor orientation_rotation {movetype_constructor<OrientationRotationMCMovetype>};
        MovetypeConstructor staple_exchange {movetype_constructor<StapleExchangeMCMovetype>};
        MovetypeConstructor cb_staple_regrowth {movetype_constructor<CBStapleRegrowthMCMovetype>};
 //       MovetypeConstructor ctcb_scaffold_regrowth {movetype_constructor<MCMovetype>};
    };

    const vector<MovetypeConstructor> movetype {
        movetype_constructor<IdentityMCMovetype>,
        movetype_constructor<OrientationRotationMCMovetype>,
        movetype_constructor<StapleExchangeMCMovetype>,
        movetype_constructor<CBStapleRegrowthMCMovetype>};
 //       MovetypeConstructor ctcb_scaffold_regrowth {movetype_constructor<MCMovetype>};
}

#endif // MOVETYPES_H
