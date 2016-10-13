//origami_system.h

#ifndef ORIGAMI_SYSTEM_H
#define ORIGAMI_SYSTEM_H

#include <vector>
#include <map>
#include <unordered_map>
#include <utility>
#include <string>
#include <valarray>

#include "utility.h"
#include "nearest_neighbour.h"
#include "hash.h"
#include "domain.h"

using std::vector;
using std::pair;
using std::unordered_map;
using std::string;
using std::valarray;

using namespace Utility;
using namespace DomainContainer;

namespace Origami{

    // For passing information between file objects and origami system
    struct Chain {
        bool operator==(Chain chain_2);
        int index;
        int identity;
        vector<VectorThree> positions;
        vector<VectorThree> orientations;
    };
    using Chains = vector<Chain>;

    class OrigamiSystem {
        // Cubic lattice domain-level resolution model of DNA origami
        public:
    
            OrigamiSystem(
                    const vector<vector<int>>& identities,
                    const vector<vector<string>>& sequences,
                    const Chains& chains,
                    double temp,
                    double lattice_site_volume,
                    double cation_M,
                    double staple_M,
                    bool cyclic);
            ~OrigamiSystem();

            // THESE NEED TO BE IMPLEMENTED TO DEAL WITH THE DOMAIN POINTER VECTOR
            OrigamiSystem(const OrigamiSystem&) = default;
            OrigamiSystem& operator=(const OrigamiSystem&) = default;
            OrigamiSystem(OrigamiSystem&&) = default;
            OrigamiSystem& operator=(OrigamiSystem&&) = default;
    
            // Configuration independent system properties
            const vector<vector<int>> m_identities;
            const vector<vector<string>> m_sequences;
            const double m_temp;
            const double m_cation_M;
            const double m_staple_M;
            const double m_volume;
            const bool m_cyclic;
            const int c_scaffold {0};
    
            // Configuration properties
            bool m_constraints_violated {false};
            vector<Domain*> get_chain(int c_i);
            inline vector<vector<Domain*>> get_chains() {return m_domains;}
            vector<Domain*> get_last_chain() {return m_domains.back();}
            Domain* get_domain(int c_i, int d_i);
            inline int num_staples() const {return m_domains.size() - 1;}
            int num_unique_staples() const;
            inline int num_domains() {return m_num_domains;};

            // CONSIDER MAKING THIS MORE ROBUST
            inline int num_bound_domain_pairs() const {return ((m_num_domains / 2) -
                    m_pos_to_unbound_d.size() / 2);}
            inline int num_fully_bound_domain_pairs() const {return m_num_fully_bound_domain_pairs;}
            inline int num_misbound_domain_pairs() const {
                    return num_bound_domain_pairs() - num_fully_bound_domain_pairs();}
            inline int num_staples_of_ident(int staple_ident) const {return
                    m_identity_to_index[staple_ident].size();}
            inline vector<int> staples_of_ident(int c_ident) {return m_identity_to_index[c_ident];}
            inline vector<int> complimentary_scaffold_domains(int staple_ident)
                    const {return m_staple_ident_to_scaffold_ds[staple_ident];}
            Chains chains() const;
            Occupancy position_occupancy(VectorThree pos) const;
            inline Domain* unbound_domain_at(VectorThree pos) const {return
                    m_pos_to_unbound_d.at(pos);};
            bool check_domains_complementary(Domain& cd_i, Domain& cd_j);
            inline double energy() const {return m_energy;}
    
            // Constraint checkers
            void check_all_constraints();
            double check_domain_constraints(
                    Domain& cd_i,
                    VectorThree pos,
                    VectorThree ore);
    
            // Configuration modifiers
            double unassign_domain(Domain& cd_i);
            int add_chain(int c_i_ident);
            int add_chain(int c_i_ident, int uc_i);
            void delete_chain(int c_i);
            double set_checked_domain_config(
                    Domain& cd_i,
                    VectorThree pos,
                    VectorThree ore);
            double set_domain_config(
                    Domain& cd_i,
                    VectorThree position,
                    VectorThree orientation);
            void set_domain_orientation(Domain& cd_i, VectorThree ore);
            void centre();

            // The index that should be assigned to the next added chain
            int m_current_c_i {};

        protected:
            virtual double bind_noncomplementary_domains(Domain& cd_i, Domain& cd_j);

        private:
            // Data
            vector<vector<Domain*>> m_domains {};
            int m_num_domains {0};

            // Keeps track of all scaffold domains complementary to a domain on
            // a given staple. Only tracks staple identity to the scaffold domain
            // indices
            vector<vector<int>> m_staple_ident_to_scaffold_ds {};

            // Position in domains array to chain index
            vector<int> m_chain_indices {};
            
            // Identity to list of all chains of that type
            vector<vector<int>> m_identity_to_index {};

            // May need to access the chain type by index in m_domains only
            vector<int> m_chain_identities {};

            // Keeps track of unbound domains but indexed by position
            unordered_map<VectorThree, Domain*> m_pos_to_unbound_d {};
            
            // The state of all positiions occupied by a domain index by position
            unordered_map<VectorThree, Occupancy> m_position_occupancies {};

            // Number of fully complimentary domains bound
            int m_num_fully_bound_domain_pairs {};

            // Energy tables index by chain/domain identity pair
            unordered_map<pair<int, int>, double> m_hybridization_energies {};
            unordered_map<pair<int, int>, double> m_stacking_energies {};

            // Current total energy of system
            double m_energy {0};
    
            // Intializers
            void initialize_complementary_associations();
            void initialize_energies();
            void initialize_config(Chains chains);

            // Accessors
            double hybridization_energy(const Domain& cd_i, const Domain& cd_j) const;
            double stacking_energy(const Domain& cd_i, const Domain& cd_j) const;
    
            // States updates
            void internal_unassign_domain(Domain& cd_i);
            double unassign_bound_domain(Domain& cd_i);
            void unassign_unbound_domain(Domain& cd_i);
            void update_domain(Domain& cd_i, VectorThree pos, VectorThree ore);
            void update_occupancies(
                    Domain& cd_i,
                    VectorThree position);
    
            // Constraint checkers
            double bind_domain(Domain& cd_i);
            double bind_complementary_domains(Domain& cd_i, Domain& cd_j);
            double check_stacking(Domain& cd_new, Domain& cd_old);
            bool check_domain_pair_constraints(Domain& cd_i);
            bool check_helical_constraints(Domain& cd_1, Domain& cd_2);
            bool check_linear_helix_rear(Domain& cd_3);
            bool check_linear_helix(VectorThree ndr_1, Domain& cd_2);
            bool check_junction_front(Domain& cd_1);
            bool check_junction_rear(Domain& cd_4);
            bool doubly_contiguous_junction(Domain& cd_1, Domain& cd_2);
            bool check_doubly_contiguous_junction(Domain& cd_2, Domain& cd_3);
            bool check_doubly_contiguous_junction(
                    Domain& cd_1,
                    Domain& cd_2,
                    Domain& cd_3,
                    Domain& cd_4);
            bool check_domain_orientations_opposing(Domain& cd_i, Domain& cd_j);
    };

    class OrigamiSystemWithoutMisbinding: public OrigamiSystem {
        public:
            using OrigamiSystem::OrigamiSystem;

        protected:
            double bind_noncomplementary_domains(Domain& cd_i, Domain& cd_j);
    };

    double molarity_to_lattice_volume(double molarity, double lattice_site_volume);
}

#endif // ORIGAMI_H
