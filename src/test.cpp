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

SCENARIO("Computed ideal random walks match hand calculated values") {
    IdealRandomWalks ideal_walk {};
    GIVEN("A starting point, endpoint, and number of remaining steps") {
        VectorThree start_pos {0, 0, 0};
        VectorThree end_pos {0, 2, 3};
        int N {5};
        double number_of_walks {10};
        THEN("The calculated value matches the hand calculated one") {
            REQUIRE(number_of_walks == ideal_walk.num_walks(start_pos, end_pos, N));
        }
        start_pos = {0, 0, 0};
        end_pos = {0, 2, 3};
        N = 6;
        number_of_walks = 0;
        THEN("The calculated value matches the hand calculated one") {
            REQUIRE(number_of_walks == ideal_walk.num_walks(start_pos, end_pos, N));
        }
        start_pos = {0, 0, 0};
        end_pos = {0, 2, 3};
        N = 7;
        number_of_walks = 665;
        THEN("The calculated value matches the hand calculated one") {
            REQUIRE(number_of_walks == ideal_walk.num_walks(start_pos, end_pos, N));
        }
        start_pos = {0, 0, 0};
        end_pos = {0, 2, 3};
        N = 51;
        number_of_walks = 5.947398897268465e+36;
        THEN("The calculated value matches the hand calculated one") {
            REQUIRE(number_of_walks == ideal_walk.num_walks(start_pos, end_pos, N));
        }
        start_pos = {1, 2, 3};
        end_pos = {3, 4, 5};
        N = 8;
        double num_walks_1 {ideal_walk.num_walks(start_pos, end_pos, N)};
        double num_walks_2 {ideal_walk.num_walks(start_pos, end_pos, N)};
        VectorThree origin {0, 0, 0};
        VectorThree DR {2, 2, 2};
        double num_walks_3 {ideal_walk.num_walks(origin, DR, N)};
        THEN("Switching start and endpoint leads to no change in num walks") {
            REQUIRE(num_walks_1 == num_walks_2);
            REQUIRE(num_walks_1 == num_walks_3);
        }
        vector<double> list_num_walks {};
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, 2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {-1, 2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, -2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, 2, -3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {-1, -2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {-1, 2, -3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, -2, -3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {-1, -2, -3}, 8));
        THEN("The sign on the value of the difference of the compenents is irrelevant") {
            for (auto n_walks: list_num_walks) {
                REQUIRE(n_walks == list_num_walks[0]);
            }
        }
        list_num_walks.clear();
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, 2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {2, 1, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {3, 2, 1}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, 3, 2}, 8));
        THEN("The order of the components of the difference is irrelevant") {
            for (auto n_walks: list_num_walks) {
                REQUIRE(n_walks == list_num_walks[0]);
            }
        }
    }
}

SCENARIO("Computed growthpoints, endpoints, and regrowth staples match expected") {
    double temp {1};
    double staple_M {1};
    double cation_M {1};
    double lattice_site_volume {1};
    bool cyclic {false};
    IdealRandomWalks ideal_random_walks {};
    GIVEN("Fully assembled snodin origami system") {
        // System setup
        string system_filename {"tests/snodin_assembled.json"};

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
        vector<Domain*> scaffold_domains {origami.get_chain(0)};
        WHEN("Calculate constraintpoints for first three domains") {
            vector<Domain*> selected_scaffold_domains {scaffold_domains.begin(),
                    scaffold_domains.begin() + 3};
            Constraintpoints constraintpoints {origami, ideal_random_walks};
            constraintpoints.calculate_constraintpoints(selected_scaffold_domains);
            THEN("Expected") {
                bool s1_is_gp {constraintpoints.is_growthpoint(selected_scaffold_domains[0])};
                REQUIRE(s1_is_gp == true);
                bool s2_is_gp {constraintpoints.is_growthpoint(selected_scaffold_domains[1])};
                REQUIRE(s2_is_gp == false);
                bool s3_is_gp {constraintpoints.is_growthpoint(selected_scaffold_domains[2])};
                REQUIRE(s3_is_gp == false);
                Domain* dtg {constraintpoints.get_domain_to_grow(selected_scaffold_domains[0])};
                int staple_c_i {origami.staples_of_ident(1)[0]};
                Domain* domain_to_grow {origami.get_chain(staple_c_i)[0]};
                REQUIRE(dtg == domain_to_grow);
                auto active_endpoints {constraintpoints.get_active_endpoints(0)};
                REQUIRE(active_endpoints.size() == 1);
                REQUIRE(active_endpoints[0].first == 2);
                Domain* sd {origami.get_domain(2, 1)};
                Domain* inactive_endpoint {constraintpoints.get_inactive_endpoints(sd)};
                REQUIRE(inactive_endpoint == selected_scaffold_domains[1]);
            }
        }
        WHEN("Calculate constraintpoints for last four domains") {
            vector<Domain*> selected_scaffold_domains {scaffold_domains.end() - 4,
                    scaffold_domains.end()};
            Constraintpoints constraintpoints {origami, ideal_random_walks};
            constraintpoints.calculate_constraintpoints(selected_scaffold_domains);
            THEN("Expected") {
                bool s1_is_gp {constraintpoints.is_growthpoint(selected_scaffold_domains[0])};
                REQUIRE(s1_is_gp == true);
                bool s2_is_gp {constraintpoints.is_growthpoint(selected_scaffold_domains[1])};
                REQUIRE(s2_is_gp == false);
                bool s3_is_gp {constraintpoints.is_growthpoint(selected_scaffold_domains[2])};
                REQUIRE(s3_is_gp == true);
                bool s4_is_gp {constraintpoints.is_growthpoint(selected_scaffold_domains[3])};
                REQUIRE(s4_is_gp == false);
                auto active_endpoints {constraintpoints.get_active_endpoints(0)};
                REQUIRE(active_endpoints.size() == 0);
                Domain* inact_d_1 {selected_scaffold_domains[1]->m_bound_domain};
                Domain* inact_d_2 {selected_scaffold_domains[3]->m_bound_domain};
                Domain* inactive_endpoint_1 {constraintpoints.get_inactive_endpoints(inact_d_1)};
                Domain* inactive_endpoint_2 {constraintpoints.get_inactive_endpoints(inact_d_2)};
                REQUIRE(inactive_endpoint_1 == selected_scaffold_domains[1]);
                REQUIRE(inactive_endpoint_2 == selected_scaffold_domains[3]);
            }
            WHEN("Constraints is updated as if growing") {
                constraintpoints.update_endpoints(selected_scaffold_domains[0]->m_bound_domain);
                constraintpoints.update_endpoints(selected_scaffold_domains[1]->m_bound_domain);
                auto active_endpoints {constraintpoints.get_active_endpoints(0)};
                REQUIRE(active_endpoints.size() == 1);
                REQUIRE(active_endpoints[0].first == selected_scaffold_domains[1]->m_d);
                constraintpoints.update_endpoints(selected_scaffold_domains[1]);
                constraintpoints.update_endpoints(selected_scaffold_domains[2]);
                active_endpoints = constraintpoints.get_active_endpoints(0);
                REQUIRE(active_endpoints.size() == 0);
                constraintpoints.update_endpoints(selected_scaffold_domains[2]->m_bound_domain);
                constraintpoints.update_endpoints(selected_scaffold_domains[3]->m_bound_domain);
                active_endpoints = constraintpoints.get_active_endpoints(0);
                REQUIRE(active_endpoints.size() == 1);
                REQUIRE(active_endpoints[0].first == selected_scaffold_domains[3]->m_d);
                constraintpoints.update_endpoints(selected_scaffold_domains[3]);
            }
        }
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

        Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
        Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
        Domain& staple_d_1 {*origami.get_domain(1, 0)};
        Domain& staple_d_2 {*origami.get_domain(1, 1)};

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
            Domain& staple2_d_1 {*origami.get_domain(2, 0)};
            origami.set_domain_config(staple_d_1, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple2_d_1, {0, 0, 0},
                    {-1, 0, 0})};
            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == staple_1_1_ene);
            }
        }
        WHEN("Staple 2 and staple 2 bind") {
            origami.add_chain(1);
            Domain& staple2_d_2 {*origami.get_domain(2, 1)};
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

        Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
        Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
        Domain& staple_d_1 {*origami.get_domain(1, 0)};
        Domain& staple_d_2 {*origami.get_domain(1, 1)};

        origami.unassign_domain(scaffold_d_1);
        origami.unassign_domain(scaffold_d_2);
        origami.unassign_domain(staple_d_1);
        origami.unassign_domain(staple_d_2);

        WHEN("The staple binds to it's domains in same helix configuration") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
            origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {-1, 0, 0});
            origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
            origami.set_domain_config(staple_d_2, {0, 1, 0}, {1, 0, 0});
            THEN("The configuration is allowed") {
                REQUIRE(origami.m_constraints_violated == false);
            }
        }
        WHEN("The staple binds to it's domains with any other configuration") {
            WHEN("Incorrect twist 1") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, -1, 0});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, 1, 0});
                THEN("The configuration is not allowed") {
                    REQUIRE(origami.m_constraints_violated == true);
                }
            }
            WHEN("Incorrect twist 2") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {1, 0, 0});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                origami.set_domain_config(staple_d_2, {0, 1, 0}, {-1, 0, 0});
                THEN("The configuration is not allowed") {
                    REQUIRE(origami.m_constraints_violated == true);
                }
            }
            WHEN("Incorrect twist 3") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 1, 0});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, -1, 0});
                THEN("The configuration is not allowed") {
                    REQUIRE(origami.m_constraints_violated == true);
                }
            }
            WHEN("Incorrect twist 4") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, -1});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, 0, 1});
                THEN("The configuration is not allowed") {
                    REQUIRE(origami.m_constraints_violated == true);
                }
            }
            WHEN("Incorrect twist 5") {
                origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                origami.set_domain_config(staple_d_1, {0, 0, 0}, {-1, 0, 0});
                origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, 0, -1});
                THEN("The configuration is not allowed") {
                    REQUIRE(origami.m_constraints_violated == true);
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

        Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
        Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
        Domain& scaffold_d_3 {*origami.get_domain(0, 2)};
        Domain& scaffold_d_4 {*origami.get_domain(0, 3)};
        Domain& staple1_d_1 {*origami.get_domain(1, 0)};
        Domain& staple1_d_2 {*origami.get_domain(1, 1)};
        Domain& staple2_d_1 {*origami.get_domain(2, 0)};
        Domain& staple2_d_2 {*origami.get_domain(2, 1)};

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
            origami.set_domain_config(staple2_d_2, {0, 0, 0}, {-1, 0, 0});
            THEN("The configuration not allowed") {
                REQUIRE(origami.m_constraints_violated == true);
            }
        }
        WHEN("The first staple is bound to it's domains as a junction") {
            origami.set_domain_config(scaffold_d_2, {0, 0, 0}, {0, 1, 0});
            origami.set_domain_config(scaffold_d_3, {0, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple2_d_1, {0, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple2_d_2, {0, 0, 0}, {0, -1, 0});
            THEN("The configuration is allowed") {
                REQUIRE(origami.m_constraints_violated == false);
            }
        }
        WHEN("The first staple is bound to it's domains with any other configuration") {
            WHEN("Same helix") {
                origami.set_domain_config(scaffold_d_2, {0, 0, 0}, {1, 0, 0});
                origami.set_domain_config(scaffold_d_3, {0, 1, 0}, {-1, 0, 0});
                origami.set_domain_config(staple2_d_1, {0, 1, 0}, {1, 0, 0});
                origami.set_domain_config(staple2_d_2, {0, 0, 0}, {-1, 0, 0});
                THEN("The configuration is allowed") {
                    REQUIRE(origami.m_constraints_violated == true);
                }
            }
        }
        WHEN("Two copies of the second staple are bound to their respective domains") {
            origami.add_chain(1);
            Domain& staple12_d_1 {*origami.get_domain(3, 0)};
            Domain& staple12_d_2 {*origami.get_domain(3, 1)};
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
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
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
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
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
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
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
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
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
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
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
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
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
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
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
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
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
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
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
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {-1, 0, 0});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
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
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
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
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
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
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
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
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
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
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
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
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
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

        for (auto domain: origami.get_chain(0)) {
            origami.unassign_domain(*domain);
        }
        origami.add_chain(1);
        origami.add_chain(2);
        WHEN("First three domains of scaffold are in same linear helix") {
            origami.set_domain_config(*origami.get_domain(0, 0), {0, 0, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.get_domain(0, 1), {0, 1, 0}, {0, 0, -1});
            origami.set_domain_config(*origami.get_domain(0, 2), {0, 2, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.get_domain(1, 0), {0, 0, 0}, {0, 0, -1});
            origami.set_domain_config(*origami.get_domain(1, 1), {0, 1, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.get_domain(2, 1), {0, 2, 0}, {0, 0, -1});
            THEN("Configuration allowed") {
                REQUIRE(origami.m_constraints_violated == false);
            }
        }
        WHEN("First three domains of scaffold are in same non-linear helix") {
            origami.set_domain_config(*origami.get_domain(0, 0), {0, 0, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.get_domain(0, 1), {0, 1, 0}, {0, 0, -1});
            origami.set_domain_config(*origami.get_domain(0, 2), {1, 1, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.get_domain(1, 0), {0, 0, 0}, {0, 0, -1});
            origami.set_domain_config(*origami.get_domain(1, 1), {0, 1, 0}, {0, 0, 1});
            origami.set_domain_config(*origami.get_domain(2, 1), {1, 1, 0}, {0, 0, -1});
            THEN("Configuration allowed") {
                REQUIRE(origami.m_constraints_violated == true);
            }
        }
    }
}
