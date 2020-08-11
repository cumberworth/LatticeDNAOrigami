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

// Shared constraint checkers
bool check_domain_orientations_opposing(Domain& cd_i, Domain& cd_j);
bool check_doubly_contig(Domain* cd_1, Domain* cd_2);
bool check_domains_exist_and_bound(vector<Domain*> cdv);
bool doubly_contiguous_helix(Domain* cd_1, Domain* cd_2);

struct DeltaConfig {
    double e {0};
    int stacked_pairs {0};
    int linear_helices {0};
    int stacked_juncts {0};
};

// Forward declaration
class OrigamiPotential;

/** Potential for fully complementary binding domains
 *
 * Includes two and three body (linear helix) terms.
 */
class BindingPotential {
  public:
    BindingPotential(OrigamiPotential& pot);
    virtual ~BindingPotential() {}

    DeltaConfig bind_domains(Domain& cd_i, Domain& cd_j);
    DeltaConfig check_stacking(Domain& cd_i, Domain& cd_j);

    bool m_constraints_violated;

  protected:
    OrigamiPotential& m_pot;
    DeltaConfig m_delta_config;

    virtual void check_constraints(Domain* cd);
    virtual bool check_regular_pair_constraints(
            Domain* cd_1,
            Domain* cd_2,
            int i);
    bool check_pair_stacked(Domain* cd_1, Domain* cd_2);
    void check_linear_helix(Domain* cd_h1, Domain* cd_h2, Domain* cd_h3);
    virtual void check_linear_helices(Domain* cd_1, Domain* cd_2, int i);
    virtual void check_central_linear_helix(Domain& cd_i, Domain& cd_j);
};

/** Adds four body terms for single and double crossover junctions */
class JunctionBindingPotential: public BindingPotential {

  public:
    using BindingPotential::BindingPotential;

  private:
    void check_constraints(Domain* cd) override;
    bool check_regular_pair_constraints(Domain* cd_1, Domain* cd_2, int i)
            override;
    void check_linear_helices(Domain* cd_1, Domain* cd_2, int i) override;
    void check_central_linear_helix(Domain& cd_i, Domain& cd_j) override;

    bool check_possible_doubly_contig_helix(Domain* cd_1, Domain* cd_2);
    void check_doubly_contig_junction(
            Domain* cd_1,
            Domain* cd_2,
            Domain* cd_3,
            Domain* cd_4);
    void check_possible_doubly_contig_junction(Domain* cd_1, Domain* cd_2);
    void check_edge_pair_junction(Domain* cd_1, Domain* cd_2, int i);

    /** Check for unstacked single junctions from first two domains
     *
     * The domains passed are the first two domains of the junction.
     * A total of nine combinations of domains will be tested,
     * including two different kink pairs
     */
    void check_forward_single_junction(Domain* cd_1, Domain* cd_2);

    /** Check for unstacked single junctions from last two domains
     *
     * The domains passed are the last two domains of the junction.
     * A total of nine combinations of domains will be tested,
     * including two different kink pairs
     */
    void check_backward_single_junction(Domain* cd_1, Domain* cd_2);

    /** Check for unstacked single junctions from kink pair
     *
     * The domains passed are the kink pair. A total of nine
     * combinations of domains will be tested.
     */
    void check_central_single_junction(Domain* cd_1, Domain* cd_2);

    /** Check for unstacked single junctions
     *
     * This will subtract a stacked pair from configuration that are
     * implied to be less stacked than the other rules would imply.
     */
    void check_single_junction(
            Domain* cd_j1,
            Domain* cd_j2,
            Domain* cd_j3,
            Domain* cd_j4,
            Domain* cd_k1,
            Domain* cd_k2);
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
