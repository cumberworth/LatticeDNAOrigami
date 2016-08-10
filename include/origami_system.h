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

using std::vector;
using std::pair;
using std::unordered_map;
using std::string;
using std::valarray;

using namespace Utility;

namespace Origami{

    struct Chain {
        int index;
        int identity;
        vector<VectorThree> positions;
        vector<VectorThree> orientations;
    };
    using Chains = vector<Chain>;

    struct ConstraintViolation {};

    class OrigamiSystem {
        // Cubic lattice domain-level resolution model of DNA origami
        public:
    
            // Configuration independent system properties
            const vector<vector<int>> m_identities;
            const vector<vector<string>> m_sequences;
            const double m_temp;
            const double m_volume;
            const double m_cation_M;
            const double m_strand_M;
            static const int c_scaffold {0};
    
            // Constructor and destructor
            OrigamiSystem(
                    const vector<vector<int>>& identities,
                    const vector<vector<string>>& sequences,
                    const Chains& chains,
                    double temp,
                    double volume,
                    double cation_M,
                    double strand_M);
            ~OrigamiSystem() = default;
    
            // Copy and move
            OrigamiSystem(const OrigamiSystem&) = default;
            OrigamiSystem& operator=(const OrigamiSystem&) = default;
            OrigamiSystem(OrigamiSystem&&) = default;
            OrigamiSystem& operator=(OrigamiSystem&&) = default;
    
            // Configuration properties
            inline unordered_map<int, int> chain_lengths() const {return m_chain_lengths;};
            inline int num_staples() const {return m_chain_lengths.size() - 1;};
            inline int num_bound_domains() const {return m_num_fully_bound_domains;};
            inline double energy() const {return m_energy;};
    
            // Staple properties
            inline int num_staples_of_ident(int staple_ident) const {return
                    m_identity_to_index[staple_ident].size();};
            inline vector<int> complimentary_scaffold_domains(int staple_ident)
                    const {return m_staple_ident_to_scaffold_ds[
                            m_chain_identities.at(staple_ident)];};
    
            // Configuration accessors
            Chains chains() const;
            inline VectorThree domain_position(CDPair cd_i) const {return
                    m_positions.at(cd_i.c)[cd_i.d];};
            inline VectorThree domain_orientation(CDPair cd_i) const {
                return m_orientations.at(cd_i.c)[cd_i.d];};
            inline Occupancy position_occupancy(VectorThree pos) const {
                    return m_position_occupancies.at(pos);};
            inline Occupancy domain_occupancy(CDPair cd_i) const {return
                    m_domain_occupancies.at(cd_i);};
            inline CDPair domain_bound_to(CDPair cd_i) const {return
                    m_bound_d_to_bound_d.at(cd_i);};
            inline CDPair unbound_domain_at(VectorThree pos) const {return
                    m_pos_to_unbound_d.at(pos);};
    
            // Constraint checkers
            void check_all_constraints() const;
            double check_domain_constraints(
                    CDPair cd_i,
                    VectorThree pos,
                    VectorThree ore);
    
            // Configuration modifiers
            double unassign_domain(CDPair cd_i);
            int add_chain(int c_i_ident);
            int add_chain(int c_i_ident, int uc_i);
            void delete_chain(int c_i);
            void set_checked_domain_config(
                    CDPair cd_i,
                    VectorThree pos,
                    VectorThree ore);
            double set_domain_config(
                    CDPair cd_i,
                    VectorThree position,
                    VectorThree orientation);
            void set_domain_orientation(CDPair cd_i, VectorThree ore);
            void centre();

        protected:
    
            // Accessors
            double hybridization_energy(CDPair cd_i, CDPair cd_j) const;
            double stacking_energy(CDPair cd_i, CDPair cd_j) const;
            virtual CDPair increment_index(CDPair cd_i, int incr) = 0;
    
        private:
    
            // Indexing and associative variables
            vector<int> m_chain_indices;
            vector<vector<int>> m_staple_ident_to_scaffold_ds;
            unordered_map<CDPair, CDPair> m_bound_d_to_bound_d;
            unordered_map<VectorThree, CDPair> m_pos_to_unbound_d;
            unordered_map<int, int> m_chain_identities;
            vector<vector<int>> m_identity_to_index;
            int m_current_c_i;
    
            // Configuration variables
            unordered_map<int, int> m_chain_lengths;
            unordered_map<int, vector<VectorThree>> m_positions;
            unordered_map<int, vector<VectorThree>> m_orientations;
            
            // Occupancy variables
            unordered_map<VectorThree, Occupancy> m_position_occupancies;
            unordered_map<CDPair, Occupancy> m_domain_occupancies;
            int m_num_fully_bound_domains;

            // Energies variables
            //vector<vector<vector<double>>> hybridization_energies;
            unordered_map<pair<CDPair, CDPair>, double> m_hybridization_energies;
            unordered_map<pair<CDPair, CDPair>, double> m_stacking_energies;
            double m_energy {0};
    
            // Intializers
            void initialize_complementary_associations();
            void initialize_energies();
            void initialize_config(Chains chains);

            // States updates
            double unassign_bound_domain(CDPair cd_i);
            void unassign_unbound_domain(CDPair cd_i);
            void update_domain(CDPair cd_i, VectorThree pos, VectorThree ore);
            void update_occupancies(
                    CDPair cd_i,
                    VectorThree position);
    
            // Constraint checkers
            double bind_domain(CDPair cd_i);
            virtual double bind_noncomplementary_domains(CDPair cd_i, CDPair cd_j) = 0;
            double bind_complementary_domains(CDPair cd_i, CDPair cd_j);
            bool check_domains_complementary(CDPair cd_i, CDPair cd_j);
            double check_stacking(CDPair cd_new, CDPair cd_old);
            void check_domain_pair_constraints(CDPair cd_i);
            void check_helical_constraints(CDPair cd_1, CDPair cd_2);
            void check_linear_helix_rear(CDPair cd_3);
            void check_linear_helix(
                    VectorThree ndr_1,
                    VectorThree pos_2,
                    VectorThree ore_2,
                    CDPair cd_2);
            void check_junction_front(CDPair cd_1);
            void check_junction_rear(CDPair cd_4);
            bool doubly_contiguous_junction(CDPair cd_1, CDPair cd_2);
            void check_doubly_contiguous_junction(CDPair cd_2, CDPair cd_3);
            void check_doubly_contiguous_junction(
                    CDPair cd_1,
                    CDPair cd_2,
                    CDPair cd_3,
                    CDPair cd_4);
            void check_domain_orientations_opposing(CDPair cd_i, CDPair cd_j);
            virtual void check_twist_constraint(
                    VectorThree ndr,
                    VectorThree ore_1,
                    VectorThree ore_2) = 0;
    };

    template<typename Origami_T>
    class OrigamiSystemWithMisbinding: public Origami_T {
        protected:
            double bind_noncomplementary_domains(CDPair cd_i, CDPair cd_j);
    };

    template<typename Origami_T>
    class OrigamiSystemWithoutMisbinding: public Origami_T {
        protected:
            double bind_noncomplementary_domains(CDPair cd_i, CDPair cd_j);
    };

    template<typename Origami_T>
    class OrigamiSystemLinearScaffold: public Origami_T {
        public:
            CDPair increment_index(CDPair cd_i, int incr);
    };

    template<typename Origami_T>
    class OrigamiSystemCyclicScaffold: public Origami_T {
        public:
            CDPair increment_index(CDPair cd_i, int incr);
    };

    template<typename Origami_T>
    class OrigamiSystemSixteen: public Origami_T {
        private:
            void check_twist_constraint(
                    VectorThree nrd,
                    VectorThree ore_1,
                    VectorThree ore_2);
    };
}

#endif // ORIGAMI_H
