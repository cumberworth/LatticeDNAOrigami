// enumerator.h

#ifndef ENUMERATOR_H
#define ENUMERATOR_H

#include "parser.h"
#include "domain.h"
#include "origami_system.h"
#include "files.h"
#include "order_params.h"
  
using namespace Parser;
using namespace DomainContainer;
using namespace Origami;
using namespace Files;
using namespace OrderParams;

namespace Enumerator {

    void print_matrix(vector<vector<long double>> matrix, string filename);

    void enumerate_main(OrigamiSystem& origami, InputParameters& params);

    class ConformationalEnumerator {
        public:
            ConformationalEnumerator(
                    OrigamiSystem& origami_system);
            void enumerate();
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

        private:
            void enumerate_domain(Domain* domain, VectorThree p_prev);
            void set_growthpoint_domains(Domain* domain, VectorThree p_new);
            void set_comp_growthpoint_domains(Domain* domain, Domain* bound_domain, VectorThree p_new);
            void set_mis_growthpoint_domains(Domain* domain, Domain* bound_domain, VectorThree p_new);
            void set_bound_domain(Domain* domain, VectorThree p_new);
            void set_unbound_domain(Domain* domain, VectorThree p_new);
            void grow_next_domain(Domain* domain, VectorThree p_new);
            void create_domains_stack();
            void create_staple_stack(Domain* domain);
            double calc_multiplier(Domain* domain, Domain* other_domain);
            int count_involved_staples(Domain* domain);
            void calc_and_save_weights();
            void add_weight_matrix_entry();
  
            OrigamiSystem& m_origami_system;
  
            // Identity to unique indices
            unordered_map<int, vector<int>> m_identity_to_indices {};

            // Domain stack for growing
            vector<Domain*> m_domains;

            // Previous growthpoint position
            VectorThree m_prev_growthpoint_p;

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
  
    class GrowthpointEnumerator {
        public:
            GrowthpointEnumerator(
                    ConformationalEnumerator& conformational_enumerator,
                    OrigamiSystem& origami_system);
            void enumerate();

        private:
            void recurse_or_enumerate_conf(
                    int staple_ident,
                    int d_i,
                    size_t staple_length,
                    Domain* old_domain);
            bool growthpoints_repeated();

            ConformationalEnumerator& m_conformational_enumerator;
  
            vector<pair<int, int>> m_staples {}; // identity, number copies
            vector<Domain*> m_unbound_system_domains {}; // keep track of available domains
            OrigamiSystem& m_origami_system;

            // Current set of growthpoints (each growthpoint is a chain id and a domain index)
            vector<pair<pair<int, int>, pair<int, int>>> m_growthpoints {};

            // Set of all previously enumerated set of growthpoints
            vector<vector<pair<pair<int, int>, pair<int, int>>>> m_enumerated_growthpoints {};
    };
  
    ConformationalEnumerator enumerate_two_domain_scaffold(OrigamiSystem& origami);
    ConformationalEnumerator enumerate_four_domain_scaffold(OrigamiSystem& origami);

}

#endif // ENUMERATOR_H
