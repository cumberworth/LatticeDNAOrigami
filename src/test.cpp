// main.cpp

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#define private public
#define protected public

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

template<typename T>
bool compare_approx_vector_contents(vector<T> first_vector, vector<T> second_vector) {
    // Check if elements the same regardless of order
    if (first_vector.size() != second_vector.size()) {
        return false;
    }
    for (auto ele1: first_vector) {
        bool v2_contains_ele {false};
        for (auto ele2: second_vector) {
            if (Approx(ele1) == ele2) {
                v2_contains_ele = true;
                break;
            }
        }
        if (not v2_contains_ele) {
            return false;
        }
    }
    return true;
}

template<typename T>
bool compare_vector_contents(vector<T> first_vector, vector<T> second_vector) {
    // Check if elements the same regardless of order
    if (first_vector.size() != second_vector.size()) {
        return false;
    }
    for (auto ele1: first_vector) {
        bool v2_contains_ele {false};
        for (auto ele2: second_vector) {
            if (ele1 == ele2) {
                v2_contains_ele = true;
                break;
            }
        }
        if (not v2_contains_ele) {
            return false;
        }
    }
    return true;
}

SCENARIO("Longest contiguous complementary sequence extracted") {
    GIVEN("A non-fully complementary pair of sequences of the same length") {
        string seq1 {"ATCGAAAAAAAAACTAA"};
        string seq2 {"TTAGAAAAACGATAAAA"};
        WHEN("The order of the sequences input is reversed.") {
            vector<string> comp_seqs1 {find_longest_contig_complement(seq1, seq2)};
            vector<string> expected_comp_seqs1 {"ATCG", "CTAA"};
            vector<string> comp_seqs2 {find_longest_contig_complement(seq2, seq1)};
            vector<string> expected_comp_seqs2 {"TTAG", "CGAT"};
            REQUIRE(comp_seqs1 == expected_comp_seqs1);
            REQUIRE(comp_seqs2 == expected_comp_seqs2);
            THEN("The associated set of sequences have the average hybridization energy") {
                double temp {350};
                double cation_M {1};
                double ene1 {0};
                int N {0};
                for (auto seq: comp_seqs1) {
                    ene1 += calc_hybridization_energy(seq, temp, cation_M);
                    N++;
                }
                ene1 /= N;

                double ene2 {0};
                N = 0;
                for (auto seq: comp_seqs2) {
                    ene2 += calc_hybridization_energy(seq, temp, cation_M);
                    N++;
                }
                ene2 /= N;
                REQUIRE(Approx(ene1) == ene2);
            }
        }
    }
    GIVEN("Two scaffold sequences from the snodin system that used to fail") {
        string seq1 {"CCTTTTTTTCTTTATA"};
        string seq2 {"TCGCTTCCTACTCCCA"};
        vector<string> comp_seqs1 {find_longest_contig_complement(seq1, seq2)};
        vector<string> expected_comp_seqs1 {"TA", "TA"};
        vector<string> comp_seqs2 {find_longest_contig_complement(seq2, seq1)};
        vector<string> expected_comp_seqs2 {"TA", "TA"};
        REQUIRE(comp_seqs1 == expected_comp_seqs1);
        REQUIRE(comp_seqs2 == expected_comp_seqs2);
    }
}

SCENARIO("Computed ideal random walks match hand calculated values") {
    IdealRandomWalks ideal_walk {};
    GIVEN("A starting point, endpoint, and number of remaining steps") {
        VectorThree start_pos {0, 0, 0};
        VectorThree end_pos {0, 2, 3};
        int N {5};
        long double number_of_walks {10};
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
            REQUIRE(Approx(number_of_walks) == ideal_walk.num_walks(start_pos, end_pos, N));
        }
        start_pos = {1, 2, 3};
        end_pos = {3, 4, 5};
        N = 8;
        long double num_walks_1 {ideal_walk.num_walks(start_pos, end_pos, N)};
        long double num_walks_2 {ideal_walk.num_walks(start_pos, end_pos, N)};
        VectorThree origin {0, 0, 0};
        VectorThree DR {2, 2, 2};
        long double num_walks_3 {ideal_walk.num_walks(origin, DR, N)};
        THEN("Switching start and endpoint leads to no change in num walks") {
            REQUIRE(num_walks_1 == num_walks_2);
            REQUIRE(num_walks_1 == num_walks_3);
        }
        vector<long double> list_num_walks {};
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
                REQUIRE(Approx(delta_e) == scaffold_staple_1_ene);
            }
        }
        WHEN("Scaffold 2 and staple 2 bind") {
            origami.set_domain_config(scaffold_d_2, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple_d_2, {0, 0, 0},
                    {-1, 0, 0})};
            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == scaffold_staple_2_ene);
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
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
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
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
                    }
                }
            }
            WHEN("Add scaffold 2 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
                    }
                }
            }
            WHEN("Add scaffold 3 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
                    }
                }
            }
            WHEN("Add scaffold 4 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {0, 0, -1});
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
                    }
                }
            }
            WHEN("Add staple1 1 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {0, 0, 1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
                    }
                }
            }
            WHEN("Add staple11 2 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {0, 0, 1});
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
                    }
                }
            }
            WHEN("Add staple2 1 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {0, 0, 1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 2, 1}, {0, 0, 1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    THEN("The configuration is not allowed") {
                        REQUIRE(origami.m_constraints_violated == true);
                    }
                }
            }
            WHEN("Add staple2 2 last") {
                WHEN("Correctly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 0, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {0, 0, 1}, {0, 0, 1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    THEN("The configuration is allowed") {
                        REQUIRE(origami.m_constraints_violated == false);
                    }
                }
                WHEN("Incorrectly") {
                    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {0, 2, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
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

SCENARIO("Example moves work as expected") {
    // Here I test the functions that are called from attempt move, but not
    // attempt move itself.
    double temp {300};
    double staple_M {1e-3};
    double cation_M {1};
    double lattice_site_volume {4e-28};
    bool cyclic {false};
    RandomGens random_gens {};
    IdealRandomWalks ideal_random_walks {};

    // Coordination number of lattice
    double k {6};
    GIVEN("Four domain loop system") {
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

/*        WHEN("CB staple exchange (insertion) move attempted") {
            double expected_new_bias {1};

            // Remove chain to get number of insertion sites correct
            origami.delete_chain(1);

            // Set initial configuration, staple 2 bound 
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, -1, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 1, 0});
            origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 1, 0});
            origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple2_d_1, {1, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple2_d_2, {1, 0, 0}, {0, -1, 0});

            // Setup movetype
            CBStapleExchangeMCMovetype movetype {origami, random_gens, ideal_random_walks};

            // Set first staple domain
            origami.add_chain(1, 1);
            Domain& staple1_d_1 {*origami.get_domain(1, 0)};
            Domain& staple1_d_2 {*origami.get_domain(1, 1)};

            double delta_e {0};
            delta_e += origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 1, 0});

            // This is done in a method that has random selections happening
            movetype.m_bias *= k * exp(-delta_e);
            expected_new_bias *= k * exp(-delta_e);

            // Set second staple domain
            vector<pair<VectorThree, VectorThree>> configs {};
            vector<double> bfactors {};
            VectorThree p_prev {0, 0, 0};
            movetype.calc_biases(staple1_d_2, p_prev, configs, bfactors);
            vector<Domain*> domains {&staple1_d_1, &staple1_d_2};
            vector<double> weights {movetype.calc_bias(bfactors, &staple1_d_2,
                    configs, p_prev, domains)};

            delta_e = origami.hybridization_energy(scaffold_d_4, staple1_d_2);
            double Ri = exp(-delta_e) + 6*4;
            expected_new_bias *= Ri;
            vector<double> expected_weights {exp(-delta_e) / Ri, 6/Ri, 6/Ri, 6/Ri, 6/Ri};
            REQUIRE(compare_approx_vector_contents(weights, expected_weights));

            // Calculate acceptance ratio
            movetype.m_modifier = 1./2;
            double calc_ratio {movetype.calc_staple_insertion_acc_ratio(staple1_d_1.m_c_ident)};
            double expected_ratio {expected_new_bias / pow(6, 2) * pow(6, 3)};
            REQUIRE(Approx(calc_ratio) == expected_ratio);

            double calc_p_accept {movetype.m_modifier};

            // This value accounts for the staple being added before the movetype
            // was created, meaning there are two extra domains
            double expected_p_accept {6 / origami.m_volume};
            REQUIRE(calc_p_accept == expected_p_accept);
        }
        WHEN("CB staple exchange (deletion) move attempted") {

            // Set initial configuration, staple 2 bound 
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, -1, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 1, 0});
            origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 1, 0});
            origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 1, 0});
            origami.set_domain_config(staple1_d_2, {0, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple2_d_1, {1, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple2_d_2, {1, 0, 0}, {0, -1, 0});

            // Setup movetype
            CBStapleExchangeMCMovetype movetype {origami, random_gens, ideal_random_walks};

            // Regrow staple in old configuration
            movetype.m_regrow_old = true;
            vector<Domain*> staple {&staple1_d_1, &staple1_d_2};
            movetype.unassign_for_regrowth(staple);
            double delta_e {movetype.set_old_growth_point(staple1_d_1, scaffold_d_1)};
            movetype.m_bias *= k * exp(-delta_e);
            movetype.grow_staple(staple1_d_1.m_d, staple);

            double expected_new_bias {k * exp(-delta_e)};
            delta_e = origami.hybridization_energy(scaffold_d_4, staple1_d_2);
            double Ri = exp(-delta_e) + 6*4;
            expected_new_bias *= Ri;

            // Expects inverted modifier
            movetype.m_modifier = 1 / movetype.m_modifier;
            double calc_ratio {movetype.calc_staple_deletion_acc_ratio(staple1_d_1.m_c_ident)};
            double expected_ratio {pow(6, 2) / (expected_new_bias * pow(6, 3))};
            REQUIRE(Approx(calc_ratio) == expected_ratio);
        }
*/
        WHEN("CB staple regrowth move attempted") {

            // Set initial configuration
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 1});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {1, 0, 0});
            origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 0, -1});
            origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple2_d_1, {1, 1, 0}, {0, 0, 1});
            origami.set_domain_config(staple2_d_2, {2, 1, 0}, {1, 0, 0});
            origami.set_domain_config(staple1_d_1, {2, 1, 0}, {1, 0, 0});
            origami.set_domain_config(staple1_d_2, {1, 0, 0}, {1, 0, 0});

            // Setup movetype
            CBStapleRegrowthMCMovetype movetype {origami, random_gens, ideal_random_walks};

            // Grow staple
            // Set growthpoint
            vector<Domain*> staple {&staple1_d_1, &staple1_d_2};
            movetype.unassign_domains(staple);
            double delta_e {origami.set_domain_config(staple1_d_2, {0, 0, 0}, {0, -1, 0})};
            double expected_new_bias {exp(-delta_e)};
            movetype.m_bias *= exp(-delta_e);

            // Grow next domain
            vector<pair<VectorThree, VectorThree>> configs {};
            vector<double> bfactors {};
            VectorThree p_prev {0, 0, 0};
            movetype.calc_biases(staple1_d_1, p_prev, configs, bfactors);
            movetype.calc_bias(bfactors, &staple1_d_1, configs, p_prev, staple);
            delta_e = origami.set_domain_config(staple1_d_1, {0, 1, 0}, {-1, 0, 0});
            double delta_e_alt = origami.hybridization_energy(staple1_d_1, scaffold_d_2);
            expected_new_bias *= 6*exp(-delta_e) + 6*exp(-delta_e_alt) + 4*6;
            REQUIRE(expected_new_bias == movetype.m_bias);

            // Regrow old
            movetype.setup_for_regrow_old();
            pair<Domain*, Domain*> growthpoint {&staple1_d_1, &staple2_d_2};
            movetype.unassign_domains(staple);
            movetype.set_growthpoint_and_grow_staple(growthpoint, staple);
            delta_e = origami.hybridization_energy(staple1_d_1, staple2_d_2);
            double expected_old_bias {exp(-delta_e) * 5*6};
            REQUIRE(expected_old_bias == movetype.m_bias);
        }
        WHEN("CTCB scaffold regrowth move attempted") {

            // Set initial config, fully bound
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, -1, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 1, 0});
            origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 1, 0});
            origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 1, 0});
            origami.set_domain_config(staple1_d_2, {0, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple2_d_1, {1, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple2_d_2, {1, 0, 0}, {0, -1, 0});

            // Setup movetype
            CTCBScaffoldRegrowthMCMovetype movetype {origami, random_gens, ideal_random_walks};
            vector<Domain*> scaffold_domains {&scaffold_d_2, &scaffold_d_3, &scaffold_d_4};
            movetype.m_constraintpoints.calculate_constraintpoints(scaffold_domains);
            set<int> staples {movetype.m_constraintpoints.staples_to_be_regrown()};
            movetype.unassign_domains({&scaffold_d_3, &scaffold_d_4});
            for (auto c_i: staples) {
                movetype.unassign_domains(origami.get_chain(c_i));
            }

            // Set S2D2
            vector<Domain*> staple2 {&staple2_d_1, &staple2_d_2};
            double delta_e {origami.set_domain_config(staple2_d_2, {1, 0, 0}, {0, -1, 0})};
            movetype.m_constraintpoints.update_endpoints(&staple2_d_2);
            movetype.m_bias *= exp(-delta_e);

            // Grow S2D1
            vector<pair<VectorThree, VectorThree>> configs {};
            vector<double> bfactors {};
            VectorThree p_prev {staple2_d_2.m_pos};
            movetype.calc_biases(staple2_d_1, p_prev, configs, bfactors);
            movetype.calc_bias(bfactors, &staple1_d_1, configs, p_prev, staple2);
            origami.set_domain_config(staple2_d_1, {1, 1, 0}, {0, -1, 0});
            movetype.m_constraintpoints.update_endpoints(&staple2_d_1);

            // Grow S3
            configs = {};
            bfactors = {};
            p_prev = scaffold_d_2.m_pos;
            movetype.calc_biases(scaffold_d_3, p_prev, configs, bfactors);
            movetype.calc_bias(bfactors, &scaffold_d_3, configs, p_prev, scaffold_domains);
            origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 1, 0});
            movetype.m_constraintpoints.update_endpoints(&scaffold_d_3);

            // Grow S4
            configs = {};
            bfactors = {};
            p_prev = scaffold_d_3.m_pos;
            movetype.calc_biases(scaffold_d_4, p_prev, configs, bfactors);
            movetype.calc_bias(bfactors, &scaffold_d_4, configs, p_prev, scaffold_domains);
            origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, -1, 0});
            movetype.m_constraintpoints.update_endpoints(&scaffold_d_4);

            // Calculate expected
            double delta_e_sd2s2d2 {origami.hybridization_energy(scaffold_d_2, staple2_d_2)};
            double delta_e_sd3s2d1 {origami.hybridization_energy(scaffold_d_3, staple2_d_1)};
            double delta_e_sd4s1d2 {origami.hybridization_energy(scaffold_d_4, staple1_d_2)};
            double expected_new_bias {exp(-delta_e_sd2s2d2) * 5*6 * exp(-delta_e_sd3s2d1)/2 *
                    exp(-delta_e_sd4s1d2)};
            REQUIRE(expected_new_bias == movetype.m_bias);

            // Regrow old
            movetype.setup_for_regrow_old();
            movetype.m_constraintpoints.reset_active_endpoints();
            movetype.unassign_domains({&scaffold_d_3, &scaffold_d_4});
            for (auto c_i: staples) {
                movetype.unassign_domains(origami.get_chain(c_i));
            }
            movetype.grow_staple_and_update_endpoints(scaffold_domains[0]);
            movetype.grow_chain(scaffold_domains);
            REQUIRE(expected_new_bias == movetype.m_bias);
        }
    }
}

SCENARIO("Connector staples are correctly identified") {
    double temp {300};
    double staple_M {1e-3};
    double cation_M {1};
    double lattice_site_volume {4e-28};
    bool cyclic {false};
    RandomGens random_gens {};
    IdealRandomWalks ideal_random_walks {};
    GIVEN("Four domain loop system") {
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
        origami.add_chain(1);
        origami.add_chain(1);

        Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
        Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
        Domain& scaffold_d_3 {*origami.get_domain(0, 2)};
        Domain& scaffold_d_4 {*origami.get_domain(0, 3)};
        Domain& staple11_d_1 {*origami.get_domain(1, 0)};
        Domain& staple11_d_2 {*origami.get_domain(1, 1)};
        Domain& staple12_d_1 {*origami.get_domain(2, 0)};
        Domain& staple12_d_2 {*origami.get_domain(2, 1)};
        Domain& staple13_d_1 {*origami.get_domain(3, 0)};
        Domain& staple13_d_2 {*origami.get_domain(3, 1)};
        vector<Domain*> staple11 {&staple11_d_1, &staple11_d_2};
        vector<Domain*> staple12 {&staple12_d_1, &staple12_d_2};
        vector<Domain*> staple13 {&staple13_d_1, &staple13_d_2};

        origami.unassign_domain(scaffold_d_1);
        origami.unassign_domain(scaffold_d_2);
        origami.unassign_domain(scaffold_d_3);
        origami.unassign_domain(scaffold_d_4);
        origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, -1, 0});
        origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 1, 0});
        origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 1, 0});
        origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, -1, 0});

        origami.unassign_domain(staple11_d_1);
        origami.unassign_domain(staple11_d_2);
        origami.unassign_domain(staple12_d_1);
        origami.unassign_domain(staple12_d_2);
        origami.unassign_domain(staple13_d_1);
        origami.unassign_domain(staple13_d_2);

        IdentityMCMovetype movetype {origami, random_gens, ideal_random_walks};

        WHEN("Case 1") {
            origami.set_domain_config(staple11_d_1, {0, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple11_d_2, {0, 0, 0}, {0, 1, 0});
            REQUIRE(movetype.staple_is_connector(staple11) == false);
        }
        WHEN("Case 2") {
            origami.set_domain_config(staple11_d_1, {1, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple11_d_2, {2, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple12_d_1, {2, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple12_d_2, {3, 1, 0}, {0, 1, 0});

            REQUIRE(movetype.staple_is_connector(staple11) == true);
        }
        WHEN("Case 3") {
            origami.set_domain_config(staple11_d_1, {1, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple11_d_2, {2, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple12_d_1, {2, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple12_d_2, {3, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple13_d_1, {3, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple13_d_2, {4, 1, 0}, {0, 1, 0});

            REQUIRE(movetype.staple_is_connector(staple12) == true);
        }
        WHEN("Case 4") {
            origami.set_domain_config(staple11_d_1, {1, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple11_d_2, {2, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple12_d_1, {2, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple12_d_2, {2, 0, 0}, {0, 1, 0});
            origami.set_domain_config(staple13_d_1, {2, 0, 0}, {0, 1, 0});
            origami.set_domain_config(staple13_d_2, {1, 0, 0}, {0, 1, 0});

            REQUIRE(movetype.staple_is_connector(staple12) == false);
        }
    }
}

SCENARIO("Origami system can reset after move") {
    double temp {300};
    double staple_M {1e-3};
    double cation_M {1};
    double lattice_site_volume {4e-28};
    bool cyclic {false};
    RandomGens random_gens {};
    IdealRandomWalks ideal_random_walks {};

    GIVEN("Four domain loop system") {
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
        vector<Domain*> staple11 {&staple1_d_1, &staple1_d_2};
        vector<Domain*> staple12 {&staple2_d_1, &staple2_d_2};

        origami.unassign_domain(scaffold_d_1);
        origami.unassign_domain(scaffold_d_2);
        origami.unassign_domain(scaffold_d_3);
        origami.unassign_domain(scaffold_d_4);
        origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, -1, 0});
        origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 1, 0});
        origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 1, 0});
        origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, -1, 0});
        origami.set_domain_config(staple1_d_1, {0, 1, 0}, {0, -1, 0});
        origami.set_domain_config(staple1_d_2, {0, 0, 0}, {0, 1, 0});
        origami.set_domain_config(staple2_d_1, {1, 0, 0}, {0, -1, 0});
        origami.set_domain_config(staple2_d_2, {1, 1, 0}, {0, 1, 0});

        WHEN("CB staple exchanges are carried out") {
            for (int i {0}; i != 10; i++) {
                CBStapleExchangeMCMovetype movetype {origami, random_gens, ideal_random_walks};
                Chains original_chains {origami.chains()};
                double original_energy {origami.energy()};
                bool accepted {movetype.attempt_move()};

                // The movetype resets to new configuration if accepted
                if (not accepted) {
                    movetype.reset_origami();
                    Chains new_chains {origami.chains()};
                    double new_energy {origami.energy()};
                    REQUIRE(compare_vector_contents(original_chains, new_chains));
                    REQUIRE(Approx(original_energy) == new_energy);
                }
            }
        }
        WHEN("CB staple regrowths are carried out") {
            for (int i {0}; i != 10; i++) {
                CBStapleRegrowthMCMovetype movetype {origami, random_gens, ideal_random_walks};
                Chains original_chains {origami.chains()};
                double original_energy {origami.energy()};
                bool accepted {movetype.attempt_move()};

                // The movetype resets to new configuration if accepted
                if (not accepted) {
                    movetype.reset_origami();
                    Chains new_chains {origami.chains()};
                    double new_energy {origami.energy()};
                    REQUIRE(compare_vector_contents(original_chains, new_chains));
                    REQUIRE(Approx(original_energy) == new_energy);
                }
            }
        }
        WHEN("CTCB scaffold regrowths are carried out") {
            for (int i {0}; i != 10; i++) {
                CTCBScaffoldRegrowthMCMovetype movetype {origami, random_gens, ideal_random_walks};
                Chains original_chains {origami.chains()};
                double original_energy {origami.energy()};
                bool accepted {movetype.attempt_move()};

                // The movetype resets to new configuration if accepted
                if (not accepted) {
                    movetype.reset_origami();
                    Chains new_chains {origami.chains()};
                    double new_energy {origami.energy()};
                    REQUIRE(compare_vector_contents(original_chains, new_chains));
                    REQUIRE(Approx(original_energy) == new_energy);
                }
            }
        }
    }
}

SCENARIO("Methods with random elements give correct distribution") {
    double eps {0.01};
    double temp {300};
    double staple_M {1};
    double cation_M {1};
    double lattice_site_volume {1};
    bool cyclic {false};
    RandomGens random_gens {};
    IdealRandomWalks ideal_random_walks {};

    GIVEN("Two domain one staple origami") {
        string system_filename {"tests/two_domain.json"};
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

        Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
        Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
        Domain& staple_d_1 {*origami.get_domain(1, 0)};
        Domain& staple_d_2 {*origami.get_domain(1, 1)};

        IdentityMCMovetype movetype {origami, random_gens, ideal_random_walks};

        GIVEN("Uniform distribution for two domain one staple system domain selection.") {
            unordered_map<pair<int, int>, double> exp_dist {{{0, 0}, 0.25}, {{0, 1}, 0.25}, 
                {{1, 0}, 0.25}, {{1, 1}, 0.25}};
            unordered_map<pair<int, int>, double> calc_dist {};
            for (int i {0}; i != 100000; i++) {
                Domain* domain {movetype.select_random_domain()};
                pair<int, int> key {domain->m_d, domain->m_c};
                calc_dist[key] += 1./100000;
            }
            bool dists_match = true;
            for (int d {0}; d != 2; d++) {
                for (int c {0}; c != 2; c++) {
                    pair<int, int> key {d, c};
                    double calc_p {calc_dist[key]};
                    double exp_p {exp_dist[key]};
                    if (calc_p < (exp_p - eps) or calc_p > (exp_p + eps)) {
                        dists_match = false;
                    }
                }
            }
            REQUIRE(dists_match == true);
        }
        GIVEN("Uniform distrubtion of positions around a previous") {
            VectorThree p_prev {3, 4, 1};
            double exp_p {1./6};
            vector<VectorThree> new_pos {{4, 4, 1}, {2, 4, 1}, {3, 5, 1},
                    {3, 3, 1}, {3, 4, 2}, {3, 4, 0}};
            unordered_map<VectorThree, double> calc_dist {};
            for (int i {0}; i != 100000; i++) {
                VectorThree pos {movetype.select_random_position(p_prev)};
                calc_dist[pos] += 1./100000;
            }
            bool dists_match = true;
            for (auto pos: new_pos) {
                double calc_p {calc_dist[pos]};
                if (calc_p < (exp_p - eps) or calc_p > (exp_p + eps)) {
                    dists_match = false;
                }
            }
            REQUIRE(dists_match == true);
        }
        GIVEN("T such that full binding very unlikely") {
            origami.unassign_domain(scaffold_d_1);
            origami.unassign_domain(scaffold_d_2);
            origami.unassign_domain(staple_d_1);
            origami.unassign_domain(staple_d_2);

            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, -1, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, -1, 0});
            origami.set_domain_config(staple_d_1, {0, 0, 0}, {0, 1, 0});
            origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, 1, 0});

            double exp_p {1./6};
            vector<VectorThree> new_pos {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
                    {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
            unordered_map<VectorThree, double> calc_dist {};
            int accepted_moves {0};
            for (int i {0}; i != 100000; i++) {
                CTCBScaffoldRegrowthMCMovetype ct_movetype {origami, random_gens, ideal_random_walks};
                bool accepted {ct_movetype.attempt_move()};
                if (not accepted) {
                    ct_movetype.reset_origami();
                }
                else {
                    origami.centre();
                    VectorThree pos {origami.get_domain(1, 1)->m_pos};
                    calc_dist[pos]++;
                    accepted_moves++;
                }
            }
            bool dists_match = true;
            for (auto pos: new_pos) {
                double calc_p {calc_dist[pos] / accepted_moves};
                if (calc_p < (exp_p - eps) or calc_p > (exp_p + eps)) {
                    dists_match = false;
                }
            }
            REQUIRE(dists_match == true);
        }
    }
}
