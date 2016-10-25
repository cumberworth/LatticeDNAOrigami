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
            MCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks) :
                    m_origami_system {origami_system},
                    m_random_gens {random_gens},
                    m_ideal_random_walks {ideal_random_walks} {};
            virtual ~MCMovetype() {};

            virtual bool attempt_move() = 0;
            void reset_origami();

            virtual string m_label() {return "MCMovetype";};

        protected:
            OrigamiSystem& m_origami_system;
            RandomGens& m_random_gens;
            IdealRandomWalks& m_ideal_random_walks;
            bool m_rejected {false};

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
            bool test_acceptance(double p_ratio);
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
            pair<Domain*, Domain*> select_new_growthpoint(vector<Domain*> selected_chain);
    };

    class MetMCMovetype: public RegrowthMCMovetype {
        public:
            using RegrowthMCMovetype::RegrowthMCMovetype;

            string m_label() {return "MetMCMovetype";};
        protected:
            double m_delta_e {0};

            void grow_chain(vector<Domain*> domains);
            void unassign_domains(vector<Domain*> domains);
    };

    class MetStapleExchangeMCMovetype: public MetMCMovetype {
        public:
            using MetMCMovetype::MetMCMovetype;
            bool attempt_move();

            string m_label() {return "MetStapleExchangeMCMovetype";}
        protected:

            // These can be overidden for a derived class the excludes misbinding
            int preconstrained_df {0};
            //DEBUG
            //int m_insertion_sites {m_origami_system.num_domains() -
            //        m_origami_system.num_bound_domain_pairs()};
            int m_insertion_sites {m_origami_system.num_domains()};

        private:
            bool staple_insertion_accepted(int c_i_ident);
            bool staple_deletion_accepted(int c_i_ident);
            bool insert_staple();
            bool delete_staple();
    };

    class MetStapleRegrowthMCMovetype: public MetMCMovetype {
        public:
            using MetMCMovetype::MetMCMovetype;
            bool attempt_move();

            string m_label() {return "MetStapleRegrowthMCMovetype";}
    };

    class CBMCMovetype: public RegrowthMCMovetype {
        public:
            using RegrowthMCMovetype::RegrowthMCMovetype;

            string m_label() {return "CBMCMovetype";};
        protected:
            double m_bias {1};
            double m_new_bias {1};
            double m_new_modifier {1};
            bool m_regrow_old {false};
            unordered_map<pair<int, int>, VectorThree> m_old_pos {};
            unordered_map<pair<int, int>, VectorThree> m_old_ore {};

            void calc_biases(
                    Domain& domain,
                    VectorThree p_prev,
                    vector<pair<VectorThree, VectorThree>>& configs,
                    vector<double>& bfactors);
            virtual vector<double> calc_bias(vector<double> bfactors,
                    Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree,
                    vector<Domain*>) = 0;
            void select_and_set_config(vector<Domain*> domains, int i);
            void select_and_set_new_config(
                    Domain& domain,
                    vector<double> weights,
                    vector<pair<VectorThree, VectorThree>> configs);
            void select_and_set_old_config(Domain& domain);
            double set_old_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            bool test_cb_acceptance();
            void unassign_domains(vector<Domain*>);
            void setup_for_regrow_old();
            vector<pair<Domain*, Domain*>> find_bound_domains(
                    vector<Domain*> selected_chain);
            pair<Domain*, Domain*> select_old_growthpoint(
                    vector<pair<Domain*, Domain*>> bound_domains);
    };

    class CBStapleExchangeMCMovetype: public CBMCMovetype {
        public:
            using CBMCMovetype::CBMCMovetype;
            CBStapleExchangeMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks) :
                    CBMCMovetype(
                            origami_system,
                            random_gens,
                            ideal_random_walks) {}
            bool attempt_move();

            string m_label() {return "CBStapleExchangeMCMovetype";};
        protected:

            // These can be overidden for a derived class the excludes misbinding
            int preconstrained_df {0};
            int m_insertion_sites {m_origami_system.num_domains() -
                    m_origami_system.num_bound_domain_pairs()};

        private:
            double calc_staple_insertion_acc_ratio(int c_i_ident);
            double calc_staple_deletion_acc_ratio(int c_i_ident);
            vector<double> calc_bias(vector<double> bfactors,
                    Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree,
                    vector<Domain*>);
            bool insert_staple();
            bool delete_staple();
            void unassign_and_delete_staple(int c_i,vector<Domain*> staple);
            void grow_chain(vector<Domain*> domains);
            void unassign_for_regrowth(vector<Domain*>);

    };

    class CBStapleRegrowthMCMovetype: public CBMCMovetype {
        public:
            using CBMCMovetype::CBMCMovetype;
            bool attempt_move();

            string m_label() {return "CBStapleRegrowthMCMovetype";};
        private:
            void set_growthpoint_and_grow_staple(
                    pair<Domain*, Domain*> growthpoint,
                    vector<Domain*> selected_chain);
            void grow_chain(vector<Domain*> domains);
            vector<double> calc_bias(vector<double> bfactors,
                    Domain*, vector<pair<VectorThree, VectorThree>>&, VectorThree,
                    vector<Domain*>);
    };

    class Constraintpoints {
        public:
            Constraintpoints(
                    OrigamiSystem& origami_system,
                    IdealRandomWalks& ideal_random_walks);

            void calculate_constraintpoints(vector<Domain*> scaffold_domains);
            inline set<int> staples_to_be_regrown() {return m_regrowth_staples;}
            bool is_growthpoint(Domain* domain);
            void add_active_endpoint(Domain* domain, VectorThree endpoint_pos);
            void reset_active_endpoints();
            void remove_active_endpoint(Domain* domain);
            void update_endpoints(Domain* domain);
            Domain* get_domain_to_grow(Domain* domain);
            bool endpoint_reached(Domain* domain, VectorThree pos);
            double calc_num_walks_prod(Domain* domain, VectorThree pos,
                    vector<Domain*> domains, int step_offset=0);

            // For debugging
            inline vector<pair<int, VectorThree>> get_active_endpoints(int c_i) {
                    return m_active_endpoints[c_i];}
            inline Domain* get_inactive_endpoints(Domain* domain) {
                return m_inactive_endpoints[domain];}

        private:
            void find_staples_growthpoints_endpoints(vector<Domain*> scaffold_domains);
            bool staple_already_checked(Domain* domain, set<int> checked_staples);
            void add_growthpoints(vector<pair<Domain*, Domain*>> potential_growthpoints);
            void add_regrowth_staples(set<int> participating_chains);
            void add_active_endpoints_on_scaffold(
                    set<int> participating_chains,
                    vector<pair<Domain*, Domain*>> potential_growthpoints);
            void scan_staple_topology(
                    Domain* domain,
                    set<int>& potential_growthpoints,
                    vector<pair<Domain*, Domain*>>& participating_chains,
                    vector<Domain*>& scaffold_domains,
                    bool& externally_bound);

            OrigamiSystem& m_origami_system;
            IdealRandomWalks& m_ideal_random_walks;

            // Scaffold domains to be regrown
            vector<Domain*> m_scaffold_domains {};

            // Staples to be regrown
            set<int> m_regrowth_staples;

            // Map from growthpoint domain to domain to grow
            unordered_map<Domain*, Domain*> m_growthpoints {};

            // Map from chain to it's active endpoints
            unordered_map<int, vector<pair<int, VectorThree>>> m_active_endpoints {};
            unordered_map<int, vector<pair<int, VectorThree>>> m_initial_active_endpoints {};

            // Map from chain to endpoints it will impose once grown
            unordered_map<Domain*, Domain*> m_inactive_endpoints {};

    };

    class CTCBScaffoldRegrowthMCMovetype: public CBMCMovetype {
        public:
            using CBMCMovetype::CBMCMovetype;
            CTCBScaffoldRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks) :
                    CBMCMovetype(
                            origami_system,
                            random_gens,
                            ideal_random_walks) {};
            bool attempt_move();

            string m_label() {return "CTCBScaffoldRegrowthMCMovetype";};
        private:
            Constraintpoints m_constraintpoints {m_origami_system, m_ideal_random_walks};
            vector<Domain*> select_scaffold_indices();
            void grow_chain(vector<Domain*> domains);
            void grow_staple_and_update_endpoints(Domain* growth_domain_old);
            vector<double> calc_bias(
                    vector<double> bfactors,
                    Domain* domain,
                    vector<pair<VectorThree, VectorThree>>& configs,
                    VectorThree p_prev,
                    vector<Domain*> domains);
    };

    // Movetype construction

    template<typename T>
    unique_ptr<MCMovetype> movetype_constructor(OrigamiSystem& origami_system,
            RandomGens& random_gens, IdealRandomWalks& ideal_random_walks);

    using MovetypeConstructor = unique_ptr<MCMovetype> (*)(OrigamiSystem&
            origami_system, RandomGens& random_gens, IdealRandomWalks& ideal_random_walks);

    struct Movetype {
        MovetypeConstructor identity {movetype_constructor<IdentityMCMovetype>};
        MovetypeConstructor orientation_rotation {movetype_constructor<OrientationRotationMCMovetype>};
        MovetypeConstructor met_staple_exchange {movetype_constructor<MetStapleExchangeMCMovetype>};
        MovetypeConstructor met_staple_regrowth {movetype_constructor<MetStapleExchangeMCMovetype>};
        MovetypeConstructor cb_staple_exchange {movetype_constructor<CBStapleExchangeMCMovetype>};
        MovetypeConstructor cb_staple_regrowth {movetype_constructor<CBStapleRegrowthMCMovetype>};
        MovetypeConstructor ctcb_scaffold_regrowth {movetype_constructor<CTCBScaffoldRegrowthMCMovetype>};
    };

    const vector<MovetypeConstructor> movetype {
        movetype_constructor<IdentityMCMovetype>,
        movetype_constructor<OrientationRotationMCMovetype>,
        movetype_constructor<MetStapleExchangeMCMovetype>,
        movetype_constructor<MetStapleRegrowthMCMovetype>,
        movetype_constructor<CBStapleExchangeMCMovetype>,
        movetype_constructor<CBStapleRegrowthMCMovetype>,
        movetype_constructor<CTCBScaffoldRegrowthMCMovetype>};
}

#endif // MOVETYPES_H
