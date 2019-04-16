// enumerator.h

#ifndef ENUMERATOR_H
#define ENUMERATOR_H

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>

#include "bias_functions.h"
#include "domain.h"
#include "order_params.h"
#include "origami_system.h"
#include "parser.h"
#include "utility.h"
  
namespace enumerator {

    using std::reference_wrapper;
    using std::pair;
    using std::string;
    using std::vector;
    using std::unordered_map;

    using biasFunctions::SystemBiases;
    using domainContainer::Domain;
    using orderParams::SystemOrderParams;
    using orderParams::OrderParam;
    using origami::OrigamiSystem;
    using parser::InputParameters;
    using utility::VectorThree;

    void print_matrix(vector<vector<long double>> matrix, string filename);

    void enumerate_main(
            OrigamiSystem& origami,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params);

    class OvercountCalculator {
        public:
            virtual double calc_multiplier(Domain* domain, Domain* other_domain) = 0;
    };

    class MisbindingOnlyOvercountCalculator: public OvercountCalculator {
        public:
            double calc_multiplier(Domain* staple_domain, Domain*);
    };

    class MaxTwoDomainOvercountCalculator: public OvercountCalculator  {
        public:
            double calc_multiplier(Domain* domain, Domain* other_domain);
            int count_involved_staples(Domain* domain);
    };

    class ConformationalEnumerator {
        public:
            ConformationalEnumerator(
                    OrigamiSystem& origami_system,
                    OvercountCalculator& m_overcount_calculator,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    vector<string> optags);
            virtual ~ConformationalEnumerator() {};
            virtual void enumerate();
            void add_staple(int staple);
            void remove_staple(int staple);
            vector<Domain*> add_growthpoint(
                    int new_c_ident,
                    int new_d_i,
                    Domain* old_domain);
            void remove_growthpoint(
                    Domain* old_domain);
            void normalize_weights();
            long double average_energy();
            long double average_bias();
            long double num_configs();
            void print_weights(string filename);

            unordered_map<Domain*, Domain*> m_growthpoints {};
  
            // Weights of states along number of staples, fully bound domain
            // pairs, and sum of distance pairs
            unordered_map<vector<int>, double> m_state_weights {};
            unordered_map<vector<int>, double> m_normalized_weights {};

        protected:
            virtual void unassign_domains(vector<vector<Domain*>> all_chains);
            void enumerate_domain(Domain* domain, VectorThree p_prev);
            void set_growthpoint_domains(Domain* domain, VectorThree p_new);
            void set_comp_growthpoint_domains(Domain* domain, Domain* bound_domain, VectorThree p_new);
            void set_mis_growthpoint_domains(Domain* domain, Domain* bound_domain, VectorThree p_new);
            void set_bound_domain(Domain* domain, VectorThree p_new);
            virtual void set_unbound_domain(Domain* domain, VectorThree p_new);
            virtual void grow_next_domain(Domain* domain, VectorThree p_new);
            virtual void create_domains_stack();
            void create_staple_stack(Domain* domain);
            void calc_and_save_weights();
            void add_weight_matrix_entry();
  
            OrigamiSystem& m_origami_system;
            OvercountCalculator& m_overcount_calculator;
            SystemOrderParams& m_ops;
            SystemBiases& m_biases;

            // Order parameters to record and output
            vector<string> m_optags;
            vector<reference_wrapper<OrderParam>> m_ops_to_output;
  
            // Identity to unique indices
            unordered_map<int, vector<int>> m_identity_to_indices {};

            // Domain stack for growing
            vector<Domain*> m_domains;

            // Previous growthpoint position stack
            vector<VectorThree> m_prev_growthpoint_ps;

            // Total system energies
            long double m_energy {0};
            long double m_average_energy {0};
            long double m_average_bias {0};

            // Partition function
            long double m_partition_f {0};
            long double m_num_configs {0};

            long double m_multiplier {1};
            long double m_prefix {1};

            // Number of unassigned domains indexed by domain identity
            unordered_map<int, int> m_identities_to_num_unassigned {};
    };
  
    class StapleConformationalEnumerator:
            public ConformationalEnumerator  {

        public:
            StapleConformationalEnumerator(
                    OrigamiSystem& origami_system,
                    OvercountCalculator& overcount_calculator,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    vector<string> optags);
            void enumerate();

        private:
            void unassign_domains(vector<vector<Domain*>>);
            void grow_next_domain(Domain* domain, VectorThree p_new);
            void create_domains_stack();
            void set_unbound_domain(Domain* domain, VectorThree p_new);
            void grow_off_scaffold(Domain* next_domain);

            unordered_map<Domain*, Domain*> m_inverse_growthpoints {};
    };

    class GrowthpointEnumerator {
        public:
            GrowthpointEnumerator(
                    ConformationalEnumerator& conformational_enumerator,
                    OrigamiSystem& m_origami_system);
            virtual ~GrowthpointEnumerator() {}
            void enumerate();

        protected:
            void iterate_staple_identities();
            virtual void iterate_domain_growthpoints(int staple_ident) = 0;

            ConformationalEnumerator& m_conformational_enumerator;
            OrigamiSystem& m_origami_system;
            vector<pair<int, int>> m_staples {}; // identity, number copies

            // Set of all previously enumerated set of growthpoints
            vector<vector<pair<pair<int, int>, pair<int, int>>>> m_enumerated_growthpoints {};
    };

    class MisbindingGrowthpointEnumerator: public GrowthpointEnumerator {
        public:
            MisbindingGrowthpointEnumerator(
                    ConformationalEnumerator& conformational_enumerator,
                    OrigamiSystem& origami_system);
            ~MisbindingGrowthpointEnumerator() {}
            void enumerate();

        private:
            void iterate_domain_growthpoints(int staple_ident);
            void recurse_or_enumerate_conf(
                    int staple_ident,
                    int d_i,
                    size_t staple_length,
                    Domain* old_domain);
            bool growthpoints_repeated();

            // keep track of available domains
            vector<Domain*> m_unbound_system_domains {};

            // Current set of growthpoints
            // Growthpoints are a chain id and a domain index
            vector<pair<pair<int, int>, pair<int, int>>> m_growthpoints {};
    };

    class NoMisbindingGrowthpointEnumerator: public GrowthpointEnumerator {
        public:
            NoMisbindingGrowthpointEnumerator(
                    ConformationalEnumerator& conformational_enumerator,
                    OrigamiSystem& origami_system);
            ~NoMisbindingGrowthpointEnumerator() {}
            void enumerate();

        private:
            void iterate_domain_growthpoints(int staple_ident);
            void recurse_or_enumerate_conf(
                    int staple_ident,
                    int d_i,
                    size_t staple_length,
                    Domain* old_domain);
            bool growthpoints_repeated();

             // keep track of available domains
            vector<vector<Domain*>> m_unbound_system_domains {};

            // Current set of growthpoints
            // Growthpoints are a chain id and a domain index
            vector<pair<pair<int, int>, pair<int, int>>> m_growthpoints {};
    };

    class StapleEnumerator {
        public:
            StapleEnumerator(
                    GrowthpointEnumerator& m_growthpoint_enumerator,
                    ConformationalEnumerator& conformational_enumerator,
                    OrigamiSystem& origami_system,
                    vector<int> excluded_staples);
            void enumerate(
                    int min_total_staples,
                    int max_total_staples,
                    int max_type_staples);

        private:
            GrowthpointEnumerator& m_growthpoint_enumerator;
            ConformationalEnumerator& m_conf_enumerator;
            OrigamiSystem& m_origami_system;

            vector<int> m_excluded_staples;
            int m_num_staple_types;

            int m_num_staple_combos {0};
            int m_cur_max_total_staples;
            int m_max_type_staples;

            void recurse(
                    int cur_num_staples,
                    int staple_type_i,
                    int cur_num_staples_i);

    };
}

#endif // ENUMERATOR_H
