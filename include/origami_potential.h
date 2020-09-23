// origami_potential.h

#ifndef ORIGAMI_POTENTIAL_H
#define ORIGAMI_POTENTIAL_H

#include <vector>

#include "domain.h"
#include "hash.h"
#include "nearest_neighbour.h"
#include "parser.h"
#include "utility.h"

namespace potential {

using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

using domainContainer::Domain;
using nearestNeighbour::ThermoOfHybrid;
using parser::InputParameters;
using utility::VectorThree;

bool check_domain_orientations_opposing(Domain& cd_i, Domain& cd_j);
bool check_domains_exist_and_bound(vector<Domain*> cdv);
bool doubly_contiguous(Domain* cd_1, Domain* cd_2);

// Domain order will be checked
bool check_pair_stacked(Domain* cd_1, Domain* cd_2);

// Domain order of j1 to j4 will be checked, but not k1 and k2
int check_junction_stacking_penalty(
        Domain& cd_j1,
        Domain& cd_j2,
        Domain& cd_j3,
        Domain& cd_j4,
        Domain& cd_k1,
        Domain& cd_k2);

struct DeltaConfig {
    double e {0};
    int stacked_pairs {0};
    int linear_helices {0};
    int stacked_juncts {0};
};

class OrigamiPotential;

/** Potential for fully complementary binding domains
 *
 * Intended as a parent to specific potentials. Needs to be reworked to make
 * that a reality.
 */
class BindingPotential {
  public:
    BindingPotential(OrigamiPotential& pot);
    virtual ~BindingPotential() {}

    /** Calculate full potential energy change */
    DeltaConfig bind_domains(Domain& cd_i, Domain& cd_j);

    /** Calculate stacking term energy change only */
    DeltaConfig check_stacking(Domain& cd_i, Domain& cd_j);

    bool m_constraints_violated;

  protected:
    OrigamiPotential& m_pot;
    DeltaConfig m_delta_config;

    /** Change in stacking and steric energy when changing one domain */
    virtual void calc_stacking_and_steric_terms(Domain& cd_i, Domain& cd_j) = 0;

    /** Penalty for triplets with one stacked pair
     *
     *  Order is not checked. Assumes second pair are stacked.
     */
    void check_triplet_single_stacking(
            Domain* cd_h1,
            Domain* cd_h2,
            Domain* cd_h3);

    /** Penalty for triplets with two stacked pairs
     *
     *  Order is not checked, but does not matter. Assumes both pairs are
     *  stacked.
     */
    void check_triplet_double_stacking(
            Domain* cd_h1,
            Domain* cd_h2,
            Domain* cd_h3);

    // Order is not checked, but does not matter
    void check_triply_contig_helix(Domain* cd_h1, Domain* cd_h2, Domain* cd_h3);
};

/** Potential as defined in PhD thesis of Alexander Cumberworth */
class JunctionBindingPotential: public BindingPotential {
  public:
    using BindingPotential::BindingPotential;

  private:
    void calc_stacking_and_steric_terms(Domain& cd_i, Domain& cd_j) override;

    /** Checks many stacking and steric terms that involve the given domain
     *
     *  However, not all combinations involving this domain are checked. Some
     *  that equally involve it and its bound domain are checked in the calling
     *  function. j is used to check terms involving doubly contiguous domain
     *  only once to prevent double counting.
     */
    void check_constraints(Domain* cd, int j);

    /** Checks stacking and steric terms of a pair of adjacent bound domains
     *
     *  i tells if it is looking backwards from the change domain or forwards.
     *  This determines which triplet and junction combinations to check. The
     *  order of the domains is not checked.
     */
    void check_regular_pair_constraints(Domain* cd_1, Domain* cd_2, int i);

    /** Checks stacking and steric terms adjacent doubly contig helix domains
     *
     *  i and j are as in the above functions. The order of the domains is not
     *  checked.
     */
    void check_doubly_contig_helix_pair(
            Domain* cd_1,
            Domain* cd_2,
            int i,
            int j);

    /** Checks stacking and steric terms adjacent doubly contig junction domains
     *
     *  i and j are as in the above functions. The order of the domains is not
     *  checked.
     */
    void check_doubly_contig_junction_pair(Domain* cd_1, Domain* cd_2, int j);

    /** Checks junction steric and stacking constraints
     *
     *  Domain order of k1 and k2 is checked, but not order of j1/j2 and j3/j4.
     */
    void check_junction(
            Domain* cd_j1,
            Domain* cd_j2,
            Domain* cd_j3,
            Domain* cd_j4,
            Domain* cd_k1,
            Domain* cd_k2);

    /** Check forward combinations of triplet stacking terms
     *
     *  The third domain is either the next along the chain, or either the
     *  next or the previous on the chain bound to cd_2. Assumes the domains
     *  are ordered.
     */
    void check_forward_triplet_stacking_combos(
            Domain* cd_1,
            Domain* cd_2,
            int i);

    /** Check backward combinations of triplet stacking terms
     *
     *  The third domain is either the previous along the chain, or either the
     *  next or the previous on the chain bound to cd_1. Assumes the domains
     *  are ordered.
     */
    void check_backward_triplet_stacking_combos(
            Domain* cd_1,
            Domain* cd_2,
            int i);
    /** Check triplet combinations around a bound pair
     *
     *  Does not check triplets of all three on one or the other chain, as this
     *  is check separately. There are four combinations to check.
     */
    void check_central_triplet_stacking_combos(Domain& cd_i, Domain& cd_j);

    /** Check for unstacked single junctions from last two domains
     *
     * The domains passed are the last two domains of the junction.
     * A total of nine combinations of domains will be tested,
     * including two different kink pairs. Domains assumed ordered.
     */
    void check_backward_single_junction(Domain* cd_1, Domain* cd_2);

    /** Check for unstacked single junctions from first two domains
     *
     * The domains passed are the first two domains of the junction.
     * A total of nine combinations of domains will be tested,
     * including two different kink pairs. Domains assumed ordered.
     */
    void check_forward_single_junction(Domain* cd_1, Domain* cd_2);

    /** Check for unstacked single junctions from kink pair
     *
     * The domains passed are the kink pair. A total of nine
     * combinations of domains will be tested. Domains assumed ordered.
     */
    void check_central_single_junction(Domain* cd_1, Domain* cd_2);
};

class MisbindingPotential {
  public:
    MisbindingPotential(OrigamiPotential& pot);
    virtual ~MisbindingPotential() {}

    virtual double bind_domains(Domain& cd_i, Domain& cd_j) = 0;
    bool m_constraints_violated;

  protected:
    OrigamiPotential& m_pot;
};

/**
 * Misbound domanis must have opposing orientation vectors
 */
class OpposingMisbindingPotential: public MisbindingPotential {

  public:
    using MisbindingPotential::MisbindingPotential;
    double bind_domains(Domain& cd_i, Domain& cd_j) override;
};

/**
 * Misbound domains must have opposing orientation vectors
 */
class DisallowedMisbindingPotential: public MisbindingPotential {

  public:
    using MisbindingPotential::MisbindingPotential;
    double bind_domains(Domain&, Domain&) override;
};

// Interface to origami potential
class OrigamiPotential {
  public:
    OrigamiPotential(
            const vector<vector<int>> m_identities,
            const vector<vector<string>>& sequences,
            const vector<double> enthaplies,
            const vector<double> entropies,
            InputParameters& params);
    ~OrigamiPotential();

    bool m_constraints_violated;

    void update_temp(double temp, double stacking_mult = 1);

    // Domain interactions
    DeltaConfig bind_domain(Domain& cd_i);
    bool check_domains_complementary(Domain& cd_i, Domain& cd_j);
    DeltaConfig check_stacking(Domain& cd_i, Domain& cd_j);

    // Energy calculations
    double hybridization_energy(const Domain& cd_i, const Domain& cd_j) const;
    double hybridization_enthalpy(const Domain& cd_i, const Domain& cd_j) const;
    double hybridization_entropy(const Domain& cd_i, const Domain& cd_j) const;
    double stacking_energy(const Domain& cd_i, const Domain& cd_j) const;
    double init_enthalpy() const;
    double init_entropy() const;
    double init_energy() const;

  private:
    string m_energy_filebase;
    double m_temp;
    const double m_cation_M; // Cation concentration (mol/L)
    const vector<vector<int>> m_identities; // Domain identities
    const vector<vector<string>> m_sequences; // Domain sequences
    const vector<double> m_complementary_enthalpies;
    const vector<double> m_complementary_entropies;

    // Containers for binding rules
    BindingPotential* m_binding_pot;
    MisbindingPotential* m_misbinding_pot;
    string m_stacking_pot;
    string m_hybridization_pot;
    bool m_apply_mean_field_cor;

    // Stacking energy if constant
    double m_stacking_ene {0};

    // Hybridization enthalpy and entropy if constant
    double m_binding_h;
    double m_binding_s;
    double m_misbinding_h;
    double m_misbinding_s;

    // CONSIDER DEFINING TYPE FOR THESE TABLES
    // Energy tables index by chain/domain identity pair
    unordered_map<pair<int, int>, double> m_hybridization_energies {};
    unordered_map<pair<int, int>, double> m_hybridization_enthalpies {};
    unordered_map<pair<int, int>, double> m_hybridization_entropies {};
    unordered_map<pair<int, int>, double> m_stacking_energies {};

    // Energies tables indexed by temperature
    unordered_map<double, unordered_map<pair<int, int>, double>>
            m_hybridization_energy_tables {};
    unordered_map<double, unordered_map<pair<int, int>, double>>
            m_hybridization_enthalpy_tables {};
    unordered_map<double, unordered_map<pair<int, int>, double>>
            m_hybridization_entropy_tables {};
    unordered_map<pair<double, double>, unordered_map<pair<int, int>, double>>
            m_stacking_energy_tables {};

    // Initiation enthalpy and entropy
    double m_init_enthalpy;
    double m_init_entropy;
    double m_init_energy;

    // Energy table preperation
    void get_energies();
    bool read_energies_from_file();
    void write_energies_to_file();
    void calc_energies();
    void calc_hybridization_energy(
            string seq_i,
            string seq_j,
            pair<int, int> key);
    void calc_hybridization_energy(pair<int, int> key);
    void set_hybridization_energy(pair<int, int> key);
    void calc_stacking_energy(string seq_i, string seq_j, pair<int, int> key);
    void calc_stacking_energy(pair<int, int> key);
};

} // namespace potential

#endif // ORIGAMI_POTENTIAL_H
