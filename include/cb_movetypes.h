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
            virtual bool attempt_move();
            void reset_internal();

            string m_label() {return "CTCBScaffoldRegrowthMCMovetype";};
        protected:
            Constraintpoints m_constraintpoints {m_origami_system, m_ideal_random_walks};
            int m_dir;
            vector<Domain*> select_indices(vector<Domain*>);
            pair<Domain*, Domain*> select_endpoints(vector<Domain*> domains, int min_size);
            vector<Domain*> select_cyclic_segment(Domain* start_domain,
                    Domain* end_domain);
            vector<Domain*> select_linear_segment(Domain* start_domain,
                    Domain* end_domain);
            void grow_chain(vector<Domain*> domains);
            void grow_staple_and_update_endpoints(Domain* growth_domain_old);
            vector<double> calc_bias(
                    vector<double> bfactors,
                    Domain* domain,
                    vector<pair<VectorThree, VectorThree>>& configs,
                    VectorThree p_prev,
                    vector<Domain*> domains);
    };

    class CTCBLinkerRegrowthMCMovetype: public CTCBScaffoldRegrowthMCMovetype {
        public:
            using CTCBScaffoldRegrowthMCMovetype::CTCBScaffoldRegrowthMCMovetype;
            bool attempt_move();

            string m_label() {return "CTCBLinkerRegrowthMCMovetype";};

        private:
            void reset_segment(vector<Domain*> segment, size_t last_di);
            void select_linkers(vector<Domain*> domains, vector<Domain*>& linker1, vector<Domain*>& linker2,
                    vector<Domain*>& central_segment);
            bool domains_bound_externally(vector<Domain*> domains);
            bool scan_for_external_scaffold_domain(Domain* domain,
                    vector<Domain*> domains, set<int>& participating_chains);
    };
}

#endif // CB_MOVETYPES_H
