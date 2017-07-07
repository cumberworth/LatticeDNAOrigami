// cb_movetypes.h

#ifndef CB_STAPLE_MOVETYPES_H
#define CB_STAPLE_MOVETYPES_H

#include "movetypes.h"
#include "top_constraint_points.h"

namespace CBMovetypes {

    using namespace Movetypes;
    using namespace TopConstraintPoints;

    class CBMCMovetype: public RegrowthMCMovetype {
        public:
            using RegrowthMCMovetype::RegrowthMCMovetype;

            virtual void reset_internal();
            string m_label() {return "CBMCMovetype";};
        protected:
            long double m_bias {1};
            long double m_new_bias {1};
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
            //double set_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            double set_old_growth_point(Domain& growth_domain_new, Domain& growth_domain_old);
            bool test_cb_acceptance();
            void unassign_domains(vector<Domain*>);
            void setup_for_regrow_old();
            vector<pair<Domain*, Domain*>> find_bound_domains(
                    vector<Domain*> selected_chain);
            pair<Domain*, Domain*> select_old_growthpoint(
                    vector<pair<Domain*, Domain*>> bound_domains);

            void update_bias(int sign);
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

    class CTCBScaffoldRegrowthMCMovetype: public CBMCMovetype {
        public:
            using CBMCMovetype::CBMCMovetype;
            CTCBScaffoldRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    InputParameters& params) :
                    CBMCMovetype(
                            origami_system,
                            random_gens,
                            ideal_random_walks,
                            params) {};
            bool attempt_move();
            void reset_internal();

            string m_label() {return "CTCBScaffoldRegrowthMCMovetype";};
        private:
            Constraintpoints m_constraintpoints {m_origami_system, m_ideal_random_walks};
            int m_dir;
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
}

#endif // CB_MOVETYPES_H
