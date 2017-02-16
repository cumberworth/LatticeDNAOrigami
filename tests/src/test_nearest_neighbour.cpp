#include <catch.hpp>

#define private public
#define protected public

#include <iostream>

#include "nearest_neighbour.h"

using std::cout;

using namespace NearestNeighbour;

SCENARIO("Longest contiguous complementary sequences are extracted") {

    GIVEN("A non-fully complementary pair of sequences of the same length") {
        string seq1 {"ATCGAAAAAAAAACTAA"};
        string seq2 {"TTAGAAAAACGATAAAA"};
        vector<string> comp_seqs1 {find_longest_contig_complement(seq1, seq2)};
        vector<string> comp_seqs2 {find_longest_contig_complement(seq2, seq1)};

        THEN("The extracted sets match expected") {
            vector<string> expected_comp_seqs1 {"ATCG", "CTAA"};
            vector<string> expected_comp_seqs2 {"TTAG", "CGAT"};
            REQUIRE(comp_seqs1 == expected_comp_seqs1);
            REQUIRE(comp_seqs2 == expected_comp_seqs2);
        }

        THEN("The average hybridization energy doesn't depend on which "
                "sequence is taken as the first") {
            double temp {350};
            double cation_M {1};
            double ene1 {0};
            int N {0};
            for (auto seq: comp_seqs1) {
                ene1 += calc_unitless_hybridization_energy(seq, temp, cation_M);
                N++;
            }
            ene1 /= N;

            double ene2 {0};
            N = 0;
            for (auto seq: comp_seqs2) {
                ene2 += calc_unitless_hybridization_energy(seq, temp, cation_M);
                N++;
            }
            ene2 /= N;
            REQUIRE(Approx(ene1) == ene2);
        }
    }

    // Why did they cause a fail?
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

SCENARIO("Computed free energies are compared to hand calculated values") {
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

        // Reduced free energies at 300 K and [cation] = 1 (A* = A * Beta)
        double temp {300};
        double cation_M {1};
        vector<double> energies {
                /* (((0.2) + (-7.2) + (2 * 2.2)) - 300 * ((-0.0057) +
                    (-0.0014) + (-0.0204) + (2 * 0.0069))) * 4.184 * 1000 /
                    8.3144598 / 300
                */
                2.5328725104506113,
                /* (((0.2) + (-7.2) + (2 * 2.2)) - 300 * ((-0.0057) +
                    (-0.0014) + (-0.0213) + (2 * 0.0069))) * 4.184 * 1000 /
                    8.3144598 / 300
                */
                2.9857702441073424,
                /* (((0.2) + (-8.5 - 9.8 - 8.5) + (2 * 2.2)) - 300 *
                    ((-0.0057) + (-0.0014) + (-0.0227 - 0.0244 - 0.0227) +
                    (2 * 0.0069))) * 4.184 * 1000 / 8.3144598 / 300
                */
                -5.4850947742870915,
                /* (((0.2) + (-7.8) + (1 * 2.2)) - 300 * ((-0.0057) + (-0.021) +
                    (1 * 0.0069))) * 4.184 * 1000 / 8.3144598 / 300
                */
                0.9057954673134645,
                /* (((0.2) + (-10.6) + (0 * 2.2)) - 300 * ((-0.0057) +
                    (-0.0014) + (-0.0272) + (0 * 0.0069))) * 4.184 * 1000 /
                    8.3144598 / 300
                */
                -0.18451389148978148,
                /* (((0.2) + (-7.6 - 7.2 - 7.6 - 7.2 - 8.4 - 8.5 - 7.8 - 8.4 -
                    8.2 - 7.8 - 8.5 - 8.2 - 8.4 - 10.6 - 8 - 9.8 - 8) +
                    (1 * 2.2)) - 300 * ((-0.0057) + (-21.3 - 20.4 - 21.3 -
                    21.3 - 22.4 - 22.7 - 21 -22.4 - 22.2 - 21 - 22.7 - 22.2 -
                    22.4 - 27.2 - 19.9 - 24.4 - 19.9) / 1000 + (1 * 0.0069))) *
                    4.184 * 1000 / 8.3144598 / 300
                */
                -43.193024598743925};

        THEN("The associated energies match expected") {
            for (size_t i {0}; i != seqs.size(); i++) {
                double hybrid_e {calc_unitless_hybridization_energy(seqs[i],
                        temp, cation_M)};
                REQUIRE(Approx(hybrid_e) == energies[i]);
            }
        }
    }
    GIVEN("A [cation] not equal to unity and an arbitrary pair of sequences") {

        // Consider also changing the temperature
        double temp {300};
        double cation_M {0.5};
        string seq {"AT"};
        /* (((0.2) + (-7.2) + (2 * 2.2)) - 300 * ((-0.0057) + (-0.0014) +
            (-0.0204) + (2 * 0.0069) + 0.368 * 1 * math.log(0.5) / 1000)) *
            4.184 * 1000 / 8.3144598 / 300
        */
        double energy {2.6612328678696606};
        THEN("The associated energies match expected") {
            double hybrid_e {calc_unitless_hybridization_energy(seq, temp, cation_M)};
            REQUIRE(Approx(hybrid_e) == energy);
        }
    }
}
