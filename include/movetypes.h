// movetypes.h

#ifndef MOVETYPES_H
#define MOVETYPES_H

#include <iostream>
#include <set>

#include <memory>

#include "random_gens.h"
#include "origami_system.h"
#include "ideal_random_walk.h"

using std::cout;
using std::unique_ptr;
using std::set;

using namespace Origami;
using namespace RandomGen;
using namespace IdealRandomWalk;

namespace Movetypes {

    // Movetype classes

    class MCMovetype {
        public:
            MCMovetype(OrigamiSystem& origami_system, RandomGens& random_gens) : m_origami_system
                    {origami_system}, m_random_gens {random_gens} {};
            virtual ~MCMovetype() {};

            virtual bool attempt_move() = 0;
            void reset_origami();

            virtual string m_label() {return "MCMovetype";};

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
            bool staple_is_connector(vector<Domain*> staple);
            bool scan_for_scaffold_domain(Domain*, set<int>& participating_chains);
    };

    class IdentityMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move() {return true;};

            string m_label() {return "IdentityMCMovetype";};
    };

    class OrientationRotationMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move();

            string m_label() {return "OrientationRotationMCMovetype";};
    };

    class RegrowthMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;

            string m_label() {return "RegrowthMCMovetype";};
        protected:
            double set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            void grow_staple(int d_i_index, vector<Domain*> selected_chain);
            virtual void grow_chain(vector<Domain*> domains) = 0;
    };

    class MetMCMovetype: public RegrowthMCMovetype {
        public:
            using RegrowthMCMovetype::RegrowthMCMovetype;

            string m_label() {return "MetMCMovetype";};
        protected:
            double m_delta_e {0};

            void grow_chain(vector<Domain*> domains);
    };

    class MetStapleExchangeMCMovetype: public MetMCMovetype {
        public:
            using MetMCMovetype::MetMCMovetype;
            bool attempt_move();

            string m_label() {return "MetStapleExchangeMCMovetype";};
        protected:

            // These can be overidden for a derived class the excludes misbinding
            int preconstrained_df {0};
            int m_insertion_sites {m_origami_system.num_domains()};

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

            string m_label() {return "CBMCMovetype";};
        protected:
            double m_bias {1};
            bool m_regrow_old {false};
            unordered_map<pair<int, int>, VectorThree> m_old_pos {};
            unordered_map<pair<int, int>, VectorThree> m_old_ore {};

            void calc_biases(
                    Domain& domain,
                    VectorThree p_prev,
                    vector<pair<VectorThree, VectorThree>>& configs,
                    vector<double>& bfactors);
            virtual vector<double> calc_bias(vector<double> bfactors,
                    Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree) = 0;
            void select_and_set_config(vector<Domain*> domains, int i);
            void select_and_set_new_config(
                    Domain& domain,
                    vector<double> weights,
                    vector<pair<VectorThree, VectorThree>> configs);
            void select_and_set_old_config(Domain& domain);
            double set_old_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            bool test_cb_acceptance(double new_bias);
    };

    class CBStapleExchangeMCMovetype: public CBMCMovetype {
        public:
            using CBMCMovetype::CBMCMovetype;
            bool attempt_move();

            string m_label() {return "CBStapleExchangeMCMovetype";};
        protected:

            // These can be overidden for a derived class the excludes misbinding
            int preconstrained_df {0};
            int m_insertion_sites {m_origami_system.num_domains()};

        private:
            bool staple_insertion_accepted(int c_i_ident);
            bool staple_deletion_accepted(int c_i_ident);
            vector<double> calc_bias(vector<double> bfactors,
                    Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree);
            bool insert_staple();
            bool delete_staple();
            void select_and_set_growth_point(Domain* growth_domain_new);
            void grow_chain(vector<Domain*> domains);

    };

    class CBStapleRegrowthMCMovetype: public CBMCMovetype {
        public:
            using CBMCMovetype::CBMCMovetype;
            bool attempt_move();

            string m_label() {return "CBStapleRegrowthMCMovetype";};
        private:
            void grow_chain(vector<Domain*> domains);
            vector<double> calc_bias(vector<double> bfactors,
                    Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree);
    };

    class CTCBScaffoldRegrowthMCMovetype: public CBMCMovetype {
        public:
            using CBMCMovetype::CBMCMovetype;
            bool attempt_move();

            string m_label() {return "CTCBScaffoldRegrowthMCMovetype";};
        private:
            IdealRandomWalks m_ideal_random_walks {};
            unordered_map<Domain*, Domain*> m_growthpoints {};

            // Map from chain to it's active endpoints
            unordered_map<int, vector<pair<int, VectorThree>>> m_active_endpoints {};

            // Map from chain to endpoints it will impose once grown
            unordered_map<int, vector<Domain*>> m_inactive_endpoints {};

            // Indexed by chain, then pairs of domain index and position
            vector<Domain*> select_scaffold_indices();
            set<int> find_staples_growthpoints_endpoints(vector<Domain*> scaffold_domains);
            void scan_staple_topology(
                    Domain* domain,
                    set<int>& participating_chains,
                    vector<pair<Domain*, Domain*>>& potential_growthpoints,
                    vector<Domain*>& scaffold_domains,
                    bool& externally_bound);
            void unassign_domains(vector<Domain*> domains);
            void grow_chain(vector<Domain*> domains);
            void grow_staple_and_update_endpoints(
                    vector<Domain*> domains,
                    int i);
            vector<double> calc_bias(
                    vector<double> bfactors,
                    Domain* domain,
                    vector<pair<VectorThree, VectorThree>>& configs,
                    VectorThree p_prev);
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
        MovetypeConstructor staple_exchange {movetype_constructor<MetStapleExchangeMCMovetype>};
        MovetypeConstructor cb_staple_exchange {movetype_constructor<CBStapleExchangeMCMovetype>};
        MovetypeConstructor cb_staple_regrowth {movetype_constructor<CBStapleRegrowthMCMovetype>};
        MovetypeConstructor ctcb_scaffold_regrowth {movetype_constructor<CTCBScaffoldRegrowthMCMovetype>};
    };

    const vector<MovetypeConstructor> movetype {
        movetype_constructor<IdentityMCMovetype>,
        movetype_constructor<OrientationRotationMCMovetype>,
        movetype_constructor<MetStapleExchangeMCMovetype>,
        movetype_constructor<CBStapleExchangeMCMovetype>,
        movetype_constructor<CBStapleRegrowthMCMovetype>,
        movetype_constructor<CTCBScaffoldRegrowthMCMovetype>};
}

#endif // MOVETYPES_H
