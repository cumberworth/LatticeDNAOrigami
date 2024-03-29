// origami_system.h

#ifndef ORIGAMI_SYSTEM_H
#define ORIGAMI_SYSTEM_H

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <valarray>
#include <vector>

#include "boost/serialization/vector.hpp"

#include "LatticeDNAOrigami/domain.hpp"
#include "LatticeDNAOrigami/hash.hpp"
#include "LatticeDNAOrigami/origami_potential.hpp"
#include "LatticeDNAOrigami/parser.hpp"
#include "LatticeDNAOrigami/utility.hpp"

// Forward declaration

namespace orderParams {
class SystemOrderParams;
}
namespace biasFunctions {
class SystemBiases;
}

namespace origami {

using std::string;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

using domainContainer::Domain;
using parser::InputParameters;
using potential::DeltaConfig;
using potential::OrigamiPotential;
using utility::Occupancy;
using utility::VectorThree;

// For passing information between file objects and origami system
struct Chain {
    bool operator==(Chain chain_2);
    int index;
    int identity;
    vector<VectorThree> positions;
    vector<VectorThree> orientations;

  private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version) {
        ar& index;
        ar& identity;
        ar& positions;
        ar& orientations;
    }
};
using Chains = vector<Chain>;

// Cubic lattice domain-level resolution model of DNA origami
class OrigamiSystem {
  public:
    // Standard methods
    OrigamiSystem(
            const vector<vector<int>>& identities,
            const vector<vector<string>>& sequences,
            const vector<double> enthalpies,
            const vector<double> entropies,
            const Chains& chains,
            bool cyclic,
            double staple_M,
            InputParameters& params);
    virtual ~OrigamiSystem();

    // THESE NEED TO BE IMPLEMENTED TO DEAL WITH THE DOMAIN POINTER VECTOR
    //    OrigamiSystem(const OrigamiSystem&) = default;
    //    OrigamiSystem& operator=(const OrigamiSystem&) = default;
    OrigamiSystem(OrigamiSystem&&) = default;
    //    OrigamiSystem& operator=(OrigamiSystem&&) = default;

    // Configuration independent system properties
    const vector<vector<int>> m_identities; // Domain identities
    const vector<vector<string>> m_sequences; // Domain sequences
    double m_temp; // System temperature (K)
    const double m_cation_M; // Cation concentration (mol/L)
    double m_staple_M; // Staple concentration (mol/L)
    vector<double> m_staple_us;
    vector<double> m_reduced_staple_us;
    double m_reduced_fugacity; // Also removes staple length dependence
    const bool m_cyclic; // Cyclic scaffold
    const int c_scaffold {0}; // Unique chain index of scaffold

    // Configuration properties
    orderParams::SystemOrderParams& get_system_order_params();
    biasFunctions::SystemBiases& get_system_biases();
    vector<Domain*> get_chain(int c_i);
    vector<vector<Domain*>> get_chains();
    vector<Domain*> get_last_chain();
    Domain* get_domain(int c_i, int d_i);
    vector<int> get_staple_counts();
    int num_staples() const;
    int num_unique_staples() const;
    int num_domains();
    int num_bound_domain_pairs() const;
    int num_fully_bound_domain_pairs() const;
    int num_self_bound_domain_pairs() const;
    int num_misbound_domain_pairs() const;
    int num_stacked_domain_pairs() const;
    int num_linear_helix_trips() const;
    int num_stacked_junct_quads() const;
    int num_staples_of_ident(int staple_ident) const;
    vector<int> staples_of_ident(int c_ident);
    vector<int> complementary_scaffold_domains(int staple_ident) const;
    Chains chains() const;
    Occupancy position_occupancy(VectorThree pos) const;
    Domain* unbound_domain_at(VectorThree pos) const;
    bool check_domains_complementary(Domain& cd_i, Domain& cd_j);
    double energy() const;
    void update_enthalpy_and_entropy();
    double hybridization_enthalpy();
    double hybridization_entropy();
    double stacking_energy();
    bool configuration_fully_set();
    int num_unassigned_domains();
    double init_energy();

    // Constraint checkers
    void check_all_constraints();
    virtual double check_domain_constraints(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore);
    void check_distance_constraints();
    // int check_stacking(Domain& domain);

    // Configuration modifiers
    virtual double unassign_domain(Domain& cd_i);
    int add_chain(int c_i_ident);
    int add_chain(int c_i_ident, int uc_i);
    virtual void delete_chain(int c_i);
    void temp_reduce_staples_by_one();
    void undo_reduce_staples_by_one();
    virtual double set_checked_domain_config(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore);
    virtual double set_domain_config(
            Domain& cd_i,
            VectorThree position,
            VectorThree orientation);
    virtual void set_domain_orientation(Domain& cd_i, VectorThree ore);
    void center(int centering_domain);
    void set_all_domains();
    void set_all_domains(Chains config);
    void set_config(Chains new_config);

    // System state modifiers
    void update_temp(double temp, double stacking_mult = 1);
    void update_staple_us(double temp, double staple_u_mult);
    virtual void update_bias_mult(double) {};

    // Constraints state
    bool m_constraints_violated {false};

    // The index that should be assigned to the next added chain
    int m_current_c_i {};

  protected:
    // Bookeeping stuff, could probably organize better
    vector<vector<Domain*>> m_domains {}; // Domains grouped by chain
    string m_domain_type;
    int m_num_domains {0}; // Total domains in system
    int m_num_staples {0};
    vector<vector<int>> m_staple_ident_to_scaffold_ds {}; // Staple ID to comp
                                                          // scaffold domain i
    vector<int> m_chain_indices {}; // Working to unique index
    vector<int> m_chain_identities {}; // Working index to id
    vector<vector<int>> m_identity_to_index {}; // ID to unique indices
    unordered_map<VectorThree, Domain*>
            m_pos_to_unbound_d {}; // Position to unbound domain
    unordered_map<VectorThree, Occupancy>
            m_position_occupancies {}; // State of positions
    int m_num_bound_domain_pairs {0}; // Num bound domains pairs
    int m_num_fully_bound_domain_pairs {
            0}; // Num bound fully complementary domain pairs
    int m_num_self_bound_domain_pairs {0}; // Num self-misbound domain pairs
    int m_num_unassigned_domains {0};
    int m_num_stacked_domain_pairs {0};
    int m_num_linear_helix_trips {0};
    int m_num_stacked_junct_quads {0};

    unique_ptr<orderParams::SystemOrderParams> m_ops;
    unique_ptr<biasFunctions::SystemBiases> m_biases;

    double m_energy {0};
    double m_hyb_enthalpy;
    double m_hyb_entropy;
    double m_stacking_energy;
    bool m_apply_mean_field_cor;
    OrigamiPotential m_pot;

    // Intializers
    void initialize_complementary_associations();
    void initialize_scaffold(Chain scaffold_chain);
    void initialize_staples(Chains chain);

    // States updates
    DeltaConfig internal_unassign_domain(Domain& cd_i);
    double unassign_bound_domain(Domain& cd_i);
    void unassign_unbound_domain(Domain& cd_i);
    void update_domain(Domain& cd_i, VectorThree pos, VectorThree ore);
    void update_occupancies(Domain& cd_i, VectorThree position);
    void update_energy();

    // Constraint checkers
    DeltaConfig internal_check_domain_constraints(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore);
};

class OrigamiSystemWithBias: public OrigamiSystem {
  public:
    OrigamiSystemWithBias(
            const vector<vector<int>>& identities,
            const vector<vector<string>>& sequences,
            const vector<double> enthalpies,
            const vector<double> entropies,
            const Chains& chains,
            bool cyclic,
            double staple_M,
            InputParameters& params);

    // Constraint checkers
    double check_domain_constraints(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore);

    // Configuration modifiers
    double unassign_domain(Domain& cd_i);
    // Need to make the base one virtual still
    // int add_chain(int c_i_ident);
    // int add_chain(int c_i_ident, int uc_i);
    void delete_chain(int c_i);
    double set_checked_domain_config(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore);
    double set_domain_config(
            Domain& cd_i,
            VectorThree position,
            VectorThree orientation);
    // void set_domain_orientation(Domain& cd_i, VectorThree ore);
};

// Moved from main
OrigamiSystem* setup_origami(InputParameters& params);

// Convert concentration to reduced (u/kb) chemical potential
vector<double> molarity_to_chempots(
        double molarity,
        double temp,
        vector<vector<int>> identities);

vector<double> chempots_to_reduced_chempots(
        vector<double> staple_us,
        double temp);

} // namespace origami

#endif // ORIGAMI_SYSTEM_H
