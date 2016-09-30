// main.cpp

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include "parser.h"
#include "nearest_neighbour.h"
#include "domain.h"
#include "files.h"
#include "simulation.h"

using std::cout;

using namespace Parser;
using namespace Utility;
using namespace Files;
using namespace Origami;
using namespace DomainContainer;
using namespace Simulation;
using namespace Movetypes;
using namespace NearestNeighbour;

SCENARIO("Longest contiguous complementary sequence extracted") {
    GIVEN("A non-fully complementary pair of sequences") {
        string seq1 {"AAACGCGAAA"};
        string seq2 {"CCCGCGCCCC"};
        vector<string> comp_seqs {find_longest_contig_complement(seq1, seq2)};
        REQUIRE(comp_seqs[0] == "CGCG");
    }
}

SCENARIO("Computed energies match hand calculated values") {
    double cation_M {1};
    double temp {300};
    GIVEN("A set of fully complimentary sequences and hand calculated energies") {
        vector<string> seqs {
                // Two terminal AT (simple palindrom)
                "AT",
                // Reverse complement of previous (simple palindrom)
                "TA",
                // Palindrome
                "TGCA",
                // One terminal AT
                "CT",
                // No terminal AT
                "CG",
                // Sequence with every pair in the table
                "AATTACAGTCTGACGGCC"};

        // Energies (/K)
        vector<double> energies {
                // (((0.2) + (-7.2) + (2 * 2.2)) - 300 * ((-0.0057) + (-0.0014) + (-0.0204) + (2 * 0.0069))) * 4.184 * 1000 / 8.3144598 / 300
                2.5328725104506113,
                // (((0.2) + (-7.2) + (2 * 2.2)) - 300 * ((-0.0057) + (-0.0014) + (-0.0213) + (2 * 0.0069))) * 4.184 * 1000 / 8.3144598 / 300
                2.9857702441073424,
                // (((0.2) + (-8.5 - 9.8 - 8.5) + (2 * 2.2)) - 300 * ((-0.0057) + (-0.0014) + (-0.0227 - 0.0244 - 0.0227) + (2 * 0.0069))) * 4.184 * 1000 / 8.3144598 / 300
                -5.4850947742870915,
                // (((0.2) + (-7.8) + (1 * 2.2)) - 300 * ((-0.0057) + (-0.021) + (1 * 0.0069))) * 4.184 * 1000 / 8.3144598 / 300
                0.9057954673134645,
                // (((0.2) + (-10.6) + (0 * 2.2)) - 300 * ((-0.0057) + (-0.0014) + (-0.0272) + (0 * 0.0069))) * 4.184 * 1000 / 8.3144598 / 300
                -0.18451389148978148,
                // (((0.2) + (-7.6 - 7.2 - 7.6 - 7.2 - 8.4 - 8.5 - 7.8 - 8.4 - 8.2 - 7.8 - 8.5 - 8.2 - 8.4 - 10.6 - 8 - 9.8 - 8) + (1 * 2.2)) - 300 * ((-0.0057) + (-21.3 - 20.4 - 21.3 - 21.3 - 22.4 - 22.7 - 21 -22.4 - 22.2 - 21 - 22.7 - 22.2 - 22.4 - 27.2 - 19.9 - 24.4 - 19.9) / 1000 + (1 * 0.0069))) * 4.184 * 1000 / 8.3144598 / 300
                -43.193024598743925};

        for (size_t i {0}; i != seqs.size(); i++) {
            double hybrid_e {calc_hybridization_energy(seqs[i], temp, cation_M)};
            REQUIRE(Approx(hybrid_e) == energies[i]);
        }
    }
    GIVEN("A cation [] not equal to 1 and a pair of sequences") {
        cation_M = 0.5;
        string seq {"AT"};
        // (((0.2) + (-7.2) + (2 * 2.2)) - 300 * ((-0.0057) + (-0.0014) + (-0.0204) + (2 * 0.0069) + 0.368 * 1 * math.log(0.5) / 1000)) * 4.184 * 1000 / 8.3144598 / 300
        double energy {2.6612328678696606};
        double hybrid_e {calc_hybridization_energy(seq, temp, cation_M)};
        REQUIRE(Approx(hybrid_e) == energy);
    }
    GIVEN("An origami system with 2 scaffold domains and one 2 domain staple") {

        // System setup
        double staple_M {1e-3};
        double lattice_site_volume {4e-28};
        bool cyclic {false};
        string system_filename {"tests/two_domain-one_staple.json"};

        OrigamiInputFile origami_input {system_filename};
        vector<vector<int>> identities {origami_input.m_identities};
        vector<vector<string>> sequences {origami_input.m_sequences};
        vector<Chain> configs {origami_input.m_chains};

        OrigamiSystem origami {
                identities,
                sequences,
                configs,
                temp,
                staple_M,
                cation_M,
                lattice_site_volume,
                cyclic};

        Domain& scaffold_d_1 {*origami.m_domains[0][0]};
        Domain& scaffold_d_2 {*origami.m_domains[0][1]};
        Domain& staple_d_1 {*origami.m_domains[1][0]};
        Domain& staple_d_2 {*origami.m_domains[1][1]};

        origami.unassign_domain(scaffold_d_1);
        origami.unassign_domain(scaffold_d_2);
        origami.unassign_domain(staple_d_1);
        origami.unassign_domain(staple_d_2);

        // Hand input values (sequences are both 5' to 3', but staple domains ordered 3' to 5')
        string scaffold_seq_1 {"AATTGGCAGTTTCACA"};
        string scaffold_seq_2 {"CTGCCATACACTAATA"};
        string staple_seq_1 {"TGTGAAACTGCCAATT"};
        string staple_seq_2 {"TATTAGTGTATGGCAG"};

        // Complementary energies
        double scaffold_staple_1_ene {calc_hybridization_energy(scaffold_seq_1,
                temp, cation_M)};
        double scaffold_staple_2_ene {calc_hybridization_energy(scaffold_seq_2,
                temp, cation_M)};

        // Scaffold misbound energy
        vector<string> scaffold_1_2_seqs {find_longest_contig_complement(
                scaffold_seq_1, scaffold_seq_2)};
        double ene_sum {0};
        for (auto seq: scaffold_1_2_seqs) {
            ene_sum += calc_hybridization_energy(seq, temp, cation_M);
        }
        double scaffold_1_2_ene {ene_sum / scaffold_1_2_seqs.size()};

        // Staple misbound energies
        vector<string> staple_1_1_seqs {find_longest_contig_complement(
                staple_seq_1, staple_seq_1)};
        ene_sum = 0;
        for (auto seq: staple_1_1_seqs) {
            ene_sum += calc_hybridization_energy(seq, temp, cation_M);
        }
        double staple_1_1_ene {ene_sum / staple_1_1_seqs.size()};

        vector<string> staple_2_2_seqs {find_longest_contig_complement(
                staple_seq_2, staple_seq_2)};
        ene_sum = 0;
        for (auto seq: staple_2_2_seqs) {
            ene_sum += calc_hybridization_energy(seq, temp, cation_M);
        }
        double staple_2_2_ene {ene_sum / staple_2_2_seqs.size()};

        vector<string> staple_1_2_seqs {find_longest_contig_complement(
                staple_seq_1, staple_seq_2)};
        ene_sum = 0;
        for (auto seq: staple_1_2_seqs) {
            ene_sum += calc_hybridization_energy(seq, temp, cation_M);
        }
        double staple_1_2_ene {ene_sum / staple_1_2_seqs.size()};

        WHEN("Scaffold 1 and staple 1 bind") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple_d_1, {0, 0, 0},
                    {-1, 0, 0})};
            THEN("The change in energy is consistent") {
                REQUIRE(delta_e == scaffold_staple_1_ene);
            }
        }
        WHEN("Scaffold 2 and staple 2 bind") {
            origami.set_domain_config(scaffold_d_2, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple_d_2, {0, 0, 0},
                    {-1, 0, 0})};
            THEN("The change in energy is consistent") {
                REQUIRE(delta_e == scaffold_staple_2_ene);
            }
        }
        WHEN("Scaffold 1 and scaffold 2 bind") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(scaffold_d_2, {0, 0, 0},
                    {-1, 0, 0})};
            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == scaffold_1_2_ene);
            }
        }
        WHEN("Staple 1 and staple 2 bind") {
            origami.set_domain_config(staple_d_1, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple_d_2, {0, 0, 0},
                    {-1, 0, 0})};
            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == staple_1_2_ene);
            }
        }
        WHEN("Staple 1 and staple 1 bind") {
            origami.add_chain(1);
            Domain& staple2_d_1 {*origami.m_domains[2][0]};
            origami.set_domain_config(staple_d_1, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple2_d_1, {0, 0, 0},
                    {-1, 0, 0})};
            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == staple_1_1_ene);
            }
        }
        WHEN("Staple 2 and staple 2 bind") {
            origami.add_chain(1);
            Domain& staple2_d_2 {*origami.m_domains[2][1]};
            origami.set_domain_config(staple_d_2, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple2_d_2, {0, 0, 0},
                    {-1, 0, 0})};
            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == staple_2_2_ene);
            }
        }
    }
}

SCENARIO("Implemented helical constraints consistent with intended") {
    // This is not intended to be comprehensive; the enumeration on the two
    // domain system already ensure that the total number of configuerations is equal
    // to that calculated by hand. This will mainly be to check the junction
    // constraints, but some of the two domain configurations will be checked anyways
    double temp {300};
    double cation_M {1};
    double staple_M {1e-3};
    double lattice_site_volume {4e-28};
    bool cyclic {false};
    bool constraints_violated {false};
    GIVEN("Origami system with 2 16 nt scaffold domans and 1 2 domain staple") {

        // System setup
        string system_filename {"tests/two_domain-one_staple.json"};

        OrigamiInputFile origami_input {system_filename};
        vector<vector<int>> identities {origami_input.m_identities};
        vector<vector<string>> sequences {origami_input.m_sequences};
        vector<Chain> configs {origami_input.m_chains};

        OrigamiSystem origami {
                identities,
                sequences,
                configs,
                temp,
                staple_M,
                cation_M,
                lattice_site_volume,
                cyclic};

        Domain& scaffold_d_1 {*origami.m_domains[0][0]};
        Domain& scaffold_d_2 {*origami.m_domains[0][1]};
        Domain& staple_d_1 {*origami.m_domains[1][0]};
        Domain& staple_d_2 {*origami.m_domains[1][1]};

        origami.unassign_domain(scaffold_d_1);
        origami.unassign_domain(scaffold_d_2);
        origami.unassign_domain(staple_d_1);
        origami.unassign_domain(staple_d_2);

        WHEN("The staple binds to it's domains in same helix configuration") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
            origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {-1, 0, 0});
            origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
            origami.set_domain_config(staple_d_2, {0, 1, 0}, {1, 0, 0});
            constraints_violated = false;
            THEN("The configuration is allowed") {
                REQUIRE(constraints_violated == false);
            }
        }
        WHEN("The staple binds to it's domains with any other configuration") {
            WHEN("Incorrect twist 1") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, -1, 0});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                try {
                    origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, 1, 0});
                    constraints_violated = false;
                }
                catch (ConstraintViolation) {
                    constraints_violated = true;
                }
                THEN("The configuration is not allowed") {
                    REQUIRE(constraints_violated == true);
                }
            }
            WHEN("Incorrect twist 2") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {1, 0, 0});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                try {
                    origami.set_domain_config(staple_d_2, {0, 1, 0}, {-1, 0, 0});
                    constraints_violated = false;
                }
                catch (ConstraintViolation) {
                    constraints_violated = true;
                }
                THEN("The configuration is not allowed") {
                    REQUIRE(constraints_violated == true);
                }
            }
            WHEN("Incorrect twist 3") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 1, 0});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                try {
                    origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, -1, 0});
                    constraints_violated = false;
                }
                catch (ConstraintViolation) {
                    constraints_violated = true;
                }
                THEN("The configuration is not allowed") {
                    REQUIRE(constraints_violated == true);
                }
            }
            WHEN("Incorrect twist 4") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, -1});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                try {
                    origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, 0, 1});
                    constraints_violated = false;
                }
                catch (ConstraintViolation) {
                    constraints_violated = true;
                }
                THEN("The configuration is not allowed") {
                    REQUIRE(constraints_violated == true);
                }
            }
            WHEN("Incorrect twist 5") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                try {
                    origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, 0, -1});
                    constraints_violated = false;
                }
                catch (ConstraintViolation) {
                    constraints_violated = true;
                }
                THEN("The configuration is not allowed") {
                    REQUIRE(constraints_violated == true);
                }
            }
        }
        WHEN("The staple misbinds to the domains on the scaffold") {
            THEN("Any orientation is allowed") {
            }
        }
    }
    GIVEN("Origami system with 4 16 nt scaffold domains and 2 2 domain staples") {
        // Scaffold: 1 2 3 4, staple 1: 1 4, staple 2: 3 2

        // System setup
        string system_filename {"tests/four_domain_loop.json"};

        OrigamiInputFile origami_input {system_filename};
        vector<vector<int>> identities {origami_input.m_identities};
        vector<vector<string>> sequences {origami_input.m_sequences};
        vector<Chain> configs {origami_input.m_chains};

        OrigamiSystem origami {
                identities,
                sequences,
                configs,
                temp,
                staple_M,
                cation_M,
                lattice_site_volume,
                cyclic};

        origami.add_chain(1);
        origami.add_chain(2);

        Domain& scaffold_d_1 {*origami.m_domains[0][0]};
        Domain& scaffold_d_2 {*origami.m_domains[0][1]};
        Domain& scaffold_d_3 {*origami.m_domains[0][2]};
        Domain& scaffold_d_4 {*origami.m_domains[0][3]};
        Domain& staple1_d_1 {*origami.m_domains[1][0]};
        Domain& staple1_d_2 {*origami.m_domains[1][1]};
        Domain& staple2_d_1 {*origami.m_domains[2][0]};
        Domain& staple2_d_2 {*origami.m_domains[2][1]};

        origami.unassign_domain(scaffold_d_1);
        origami.unassign_domain(scaffold_d_2);
        origami.unassign_domain(scaffold_d_3);
        origami.unassign_domain(scaffold_d_4);
        origami.unassign_domain(staple1_d_1);
        origami.unassign_domain(staple1_d_2);
        origami.unassign_domain(staple2_d_1);
        origami.unassign_domain(staple2_d_2);

        WHEN("The second staple is bound to it's domains in the same helix") {
            origami.set_domain_config(scaffold_d_2, {0, 0, 0}, {1, 0, 0});
            origami.set_domain_config(scaffold_d_3, {0, 1, 0}, {-1, 0, 0});
            origami.set_domain_config(staple2_d_1, {0, 1, 0}, {1, 0, 0});
            try {
                origami.set_domain_config(staple2_d_2, {0, 0, 0}, {-1, 0, 0});
                constraints_violated = false;
            }
            catch (ConstraintViolation) {
                constraints_violated = true;
            }
            THEN("The configuration not allowed") {
                REQUIRE(constraints_violated == true);
            }
        }
        WHEN("The first staple is bound to it's domains as a junction") {
            origami.set_domain_config(scaffold_d_2, {0, 0, 0}, {0, 1, 0});
            origami.set_domain_config(scaffold_d_3, {0, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple2_d_1, {0, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple2_d_2, {0, 0, 0}, {0, -1, 0});
            constraints_violated = false;
            THEN("The configuration is allowed") {
                REQUIRE(constraints_violated == false);
            }
        }
        WHEN("The first staple is bound to it's domains with any other configuration") {
            WHEN("Same helix") {
                origami.set_domain_config(scaffold_d_2, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_3, {0, 1, 0}, {-1, 0, 0});
                origami.set_domain_config(staple2_d_1, {0, 1, 0}, {1, 0, 0});
                try {
                    origami.set_domain_config(staple2_d_2, {0, 0, 0}, {-1, 0, 0});
                    constraints_violated = false;
                }
                catch (ConstraintViolation) {
                    constraints_violated = true;
                }
                THEN("The configuration is allowed") {
                    REQUIRE(constraints_violated == true);
                }
            }
        }
        WHEN("Two copies of the second staple are bound to their respective domains") {
            origami.add_chain(1);
            Domain& staple12_d_1 {*origami.m_domains[3][0]};
            Domain& staple12_d_2 {*origami.m_domains[3][1]};
            WHEN("Add scaffold 1 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {1, 0, 0});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {-1, 0, 0});
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    constraints_violated = false;
                    THEN("The configuration is allowed") {
                        REQUIRE(constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    try {
                        origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                        constraints_violated = false;
                    }
                    catch (ConstraintViolation) {
                        constraints_violated = true;
                    }
                    THEN("The configuration is not allowed") {
                        REQUIRE(constraints_violated == true);
                    }
                }
            }
            WHEN("Add scaffold 2 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {1, 0, 0});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {-1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    constraints_violated = false;
                    THEN("The configuration is allowed") {
                        REQUIRE(constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    try {
                        origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                        constraints_violated = false;
                    }
                    catch (ConstraintViolation) {
                        constraints_violated = true;
                    }
                    THEN("The configuration is not allowed") {
                        REQUIRE(constraints_violated == true);
                    }
                }
            }
            WHEN("Add scaffold 3 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {1, 0, 0});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {-1, 0, 0});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    constraints_violated = false;
                    THEN("The configuration is allowed") {
                        REQUIRE(constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    try {
                        origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                        constraints_violated = false;
                    }
                    catch (ConstraintViolation) {
                        constraints_violated = true;
                    }
                    THEN("The configuration is not allowed") {
                        REQUIRE(constraints_violated == true);
                    }
                }
            }
            WHEN("Add scaffold 4 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {-1, 0, 0});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {1, 0, 0});
                    constraints_violated = false;
                    THEN("The configuration is allowed") {
                        REQUIRE(constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    try {
                        origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                        constraints_violated = false;
                    }
                    catch (ConstraintViolation) {
                        constraints_violated = true;
                    }
                    THEN("The configuration is not allowed") {
                        REQUIRE(constraints_violated == true);
                    }
                }
            }
            WHEN("Add staple1 1 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {1, 0, 0});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {-1, 0, 0});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    constraints_violated = false;
                    THEN("The configuration is allowed") {
                        REQUIRE(constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    try {
                        origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                        constraints_violated = false;
                    }
                    catch (ConstraintViolation) {
                        constraints_violated = true;
                    }
                    THEN("The configuration is not allowed") {
                        REQUIRE(constraints_violated == true);
                    }
                }
            }
            WHEN("Add staple11 2 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {1, 0, 0});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {-1, 0, 0});
                    constraints_violated = false;
                    THEN("The configuration is allowed") {
                        REQUIRE(constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    try {
                        origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                        constraints_violated = false;
                    }
                    catch (ConstraintViolation) {
                        constraints_violated = true;
                    }
                    THEN("The configuration is not allowed") {
                        REQUIRE(constraints_violated == true);
                    }
                }
            }
            WHEN("Add staple2 1 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {1, 0, 0});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {-1, 0, 0});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    constraints_violated = false;
                    THEN("The configuration is allowed") {
                        REQUIRE(constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    try {
                        origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                        constraints_violated = false;
                    }
                    catch (ConstraintViolation) {
                        constraints_violated = true;
                    }
                    THEN("The configuration is not allowed") {
                        REQUIRE(constraints_violated == true);
                    }
                }
            }
            WHEN("Add staple2 2 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {1, 0, 0});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {-1, 0, 0});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    constraints_violated = false;
                    THEN("The configuration is allowed") {
                        REQUIRE(constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    try {
                        origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                        constraints_violated = false;
                    }
                    catch (ConstraintViolation) {
                        constraints_violated = true;
                    }
                    THEN("The configuration is not allowed") {
                        REQUIRE(constraints_violated == true);
                    }
                }
            }
        }
    }
    GIVEN("Snodin origami system") {

        // System setup
        string system_filename {"tests/snodin_unbound.json"};

        OrigamiInputFile origami_input {system_filename};
        vector<vector<int>> identities {origami_input.m_identities};
        vector<vector<string>> sequences {origami_input.m_sequences};
        vector<Chain> configs {origami_input.m_chains};

        OrigamiSystem origami {
                identities,
                sequences,
                configs,
                temp,
                staple_M,
                cation_M,
                lattice_site_volume,
                cyclic};

        for (auto domain: origami.m_domains[0]) {
            origami.unassign_domain(*domain);
        }
        origami.add_chain(1);
        origami.add_chain(2);
        WHEN("First three domains of scaffold are in same linear helix") {
            origami.set_domain_config(*origami.m_domains[0][0], {0, 0, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.m_domains[0][1], {0, 1, 0}, {0, 0, -1});
            origami.set_domain_config(*origami.m_domains[0][2], {0, 2, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.m_domains[1][0], {0, 0, 0}, {0, 0, -1});
            origami.set_domain_config(*origami.m_domains[1][1], {0, 1, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.m_domains[2][1], {0, 2, 0}, {0, 0, -1});
            constraints_violated = false;
            THEN("Configuration allowed") {
                REQUIRE(constraints_violated == false);
            }
        }
        WHEN("First three domains of scaffold are in same non-linear helix") {
            origami.set_domain_config(*origami.m_domains[0][0], {0, 0, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.m_domains[0][1], {0, 1, 0}, {0, 0, -1});
            origami.set_domain_config(*origami.m_domains[0][2], {1, 1, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.m_domains[1][0], {0, 0, 0}, {0, 0, -1});
            origami.set_domain_config(*origami.m_domains[1][1], {0, 1, 0}, {0, 0, 1});
            try {
                origami.set_domain_config(*origami.m_domains[2][1], {1, 1, 0}, {0, 0, -1});
                constraints_violated = false;
            }
            catch (ConstraintViolation) {
                constraints_violated = true;
            }
            THEN("Configuration allowed") {
                REQUIRE(constraints_violated == true);
            }
        }
    }
}
