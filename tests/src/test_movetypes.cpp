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
    InputParameters params {};

    // Coordination number of lattice
    double k {6};
    GIVEN("Four domain loop system") {
        // Scaffold: 1 2 3 4, staple 1: 1 4, staple 2: 3 2

        // System setup
        string system_filename {"data/four_domain_loop.json"};

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
        SystemBias system_bias {params, origami};

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

        WHEN("CB staple regrowth move attempted") {

            // Set initial configuration
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 1});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {1, 0, 0});
            origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 0, -1});
            origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple2_d_1, {1, 1, 0}, {0, 0, 1});
            origami.set_domain_config(staple2_d_2, {2, 1, 0}, {1, 0, 0});
            origami.set_domain_config(staple1_d_1, {2, 1, 0}, {-1, 0, 0});
            origami.set_domain_config(staple1_d_2, {2, 0, 0}, {1, 0, 0});

            // Setup movetype
            CBStapleRegrowthMCMovetype movetype {origami, system_bias, random_gens,
                    ideal_random_walks, params};

            // Grow staple
            // Set growthpoint
            vector<Domain*> staple {&staple1_d_1, &staple1_d_2};
            movetype.unassign_domains(staple);
            double delta_e {origami.set_domain_config(staple1_d_2, {0, 0, 0}, {0, 0, -1})};
            double expected_new_bias {exp(-delta_e)};
            movetype.m_bias *= exp(-delta_e);

            // Grow next domain
            vector<pair<VectorThree, VectorThree>> configs {};
            vector<double> bfactors {};
            VectorThree p_prev {0, 0, 0};
            movetype.calc_biases(staple1_d_1, p_prev, configs, bfactors);
            movetype.calc_bias(bfactors, &staple1_d_1, configs, p_prev, staple);
            delta_e = origami.set_domain_config(staple1_d_1, {0, 1, 0}, {0, -1, 0});
            double delta_e_alt = origami.hybridization_energy(staple1_d_1, scaffold_d_2);
            expected_new_bias *= exp(-delta_e) + exp(-delta_e_alt) + 4*6;
            REQUIRE(Approx(expected_new_bias) == movetype.m_bias);

            // Regrow old
            movetype.setup_for_regrow_old();
            pair<Domain*, Domain*> growthpoint {&staple1_d_1, &staple2_d_2};
            movetype.unassign_domains(staple);
            movetype.set_growthpoint_and_grow_staple(growthpoint, staple);
            delta_e = origami.hybridization_energy(staple1_d_1, staple2_d_2);
            double expected_old_bias {exp(-delta_e) * (5*6)};
            REQUIRE(Approx(expected_old_bias) == movetype.m_bias);
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
            CTCBScaffoldRegrowthMCMovetype movetype {origami, system_bias, random_gens,
                    ideal_random_walks, params};
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
            long double expected_new_bias {exp(-delta_e_sd2s2d2) * 5*6 * exp(-delta_e_sd3s2d1)/2 *
                    exp(-delta_e_sd4s1d2)};
            REQUIRE(Approx(expected_new_bias) == movetype.m_bias);

            // Regrow old
            movetype.setup_for_regrow_old();
            movetype.m_constraintpoints.reset_active_endpoints();
            movetype.unassign_domains({&scaffold_d_3, &scaffold_d_4});
            for (auto c_i: staples) {
                movetype.unassign_domains(origami.get_chain(c_i));
            }
            movetype.grow_staple_and_update_endpoints(scaffold_domains[0]);
            movetype.grow_chain(scaffold_domains);
            REQUIRE(Approx(expected_new_bias) == movetype.m_bias);
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
    InputParameters params {};
    GIVEN("Four domain loop system") {
        // Scaffold: 1 2 3 4, staple 1: 1 4, staple 2: 3 2

        // System setup
        string system_filename {"data/four_domain_loop.json"};

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
        SystemBias system_bias {params, origami};

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

        IdentityMCMovetype movetype {origami, system_bias, random_gens, ideal_random_walks,
                params};

        WHEN("Case 1") {
            origami.set_domain_config(staple11_d_1, {0, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple11_d_2, {0, 0, 0}, {0, 1, 0});
            REQUIRE(movetype.staple_is_connector(staple11) == false);
        }
        WHEN("Case 2") {
            origami.set_domain_config(staple11_d_1, {1, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple11_d_2, {2, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple12_d_1, {2, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple12_d_2, {3, 1, 0}, {0, 1, 0});

            REQUIRE(movetype.staple_is_connector(staple11) == true);
        }
        WHEN("Case 3") {
            origami.set_domain_config(staple11_d_1, {1, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple11_d_2, {2, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple12_d_1, {2, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple12_d_2, {3, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple13_d_1, {3, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple13_d_2, {4, 1, 0}, {0, 1, 0});

            REQUIRE(movetype.staple_is_connector(staple12) == true);
        }
        WHEN("Case 4") {
            origami.set_domain_config(staple11_d_1, {1, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple11_d_2, {2, 1, 0}, {0, 1, 0});
            origami.set_domain_config(staple12_d_1, {2, 1, 0}, {0, -1, 0});
            origami.set_domain_config(staple12_d_2, {2, 0, 0}, {0, 1, 0});
            origami.set_domain_config(staple13_d_1, {2, 0, 0}, {0, -1, 0});
            origami.set_domain_config(staple13_d_2, {1, 0, 0}, {0, -1, 0});

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
    InputParameters params {};

    GIVEN("Four domain loop system") {
        // Scaffold: 1 2 3 4, staple 1: 1 4, staple 2: 3 2

        // System setup
        string system_filename {"data/four_domain_loop.json"};

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
        SystemBias system_bias {params, origami};

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
        origami.set_domain_config(staple1_d_1, {0, 1, 0}, {0, 1, 0});
        origami.set_domain_config(staple1_d_2, {0, 0, 0}, {0, 1, 0});
        origami.set_domain_config(staple2_d_1, {1, 0, 0}, {0, -1, 0});
        origami.set_domain_config(staple2_d_2, {1, 1, 0}, {0, -1, 0});

        WHEN("Met staple exchanges are carried out") {
            for (int i {0}; i != 10; i++) {
                MetStapleExchangeMCMovetype movetype {origami, system_bias, random_gens, ideal_random_walks,
                        params};
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
                CBStapleRegrowthMCMovetype movetype {origami, system_bias, random_gens, ideal_random_walks,
                        params};
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
                CTCBScaffoldRegrowthMCMovetype movetype {origami, system_bias, random_gens,
                        ideal_random_walks, params};
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
