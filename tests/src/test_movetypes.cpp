#include <catch.hpp>

#define private public
#define protected public

#include "parser.h"
#include "nearest_neighbour.h"
#include "domain.h"
#include "files.h"
#include "simulation.h"
#include "test_origami_system.h"

using std::cout;

using namespace Parser;
using namespace Utility;
using namespace Files;
using namespace Origami;
using namespace DomainContainer;
using namespace Simulation;
using namespace Movetypes;
using namespace NearestNeighbour;
using namespace TestOrigamiSystem;

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

// Need to do example moves for the rest of the movetypes
SCENARIO("Example moves work as expected") {
    // Here I test the functions that are called from attempt move, but not
    // attempt move itself.
    double temp {300};
    double cation_M {1};
    RandomGens random_gens {};
    IdealRandomWalks ideal_random_walks {};
    InputParameters params {};
    double k {6}; // Coordiation number of lattice

    // Tests performed on a four domain scaffold origami system
    OrigamiSystem origami {setup_four_domain_scaffold_origami(temp, cation_M)};
    origami.add_chain(1);
    origami.add_chain(2);

    // Easy reference
    Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
    Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
    Domain& scaffold_d_3 {*origami.get_domain(0, 2)};
    Domain& scaffold_d_4 {*origami.get_domain(0, 3)};
    Domain& staple1_d_1 {*origami.get_domain(1, 0)};
    Domain& staple1_d_2 {*origami.get_domain(1, 1)};
    Domain& staple2_d_1 {*origami.get_domain(2, 0)};
    Domain& staple2_d_2 {*origami.get_domain(2, 1)};

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
        CBStapleRegrowthMCMovetype movetype {origami, random_gens,
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

        THEN("Trial bias matches hand calculated value") {
            REQUIRE(Approx(expected_new_bias) == movetype.m_bias);
        }

        // Regrow old
        movetype.setup_for_regrow_old();
        pair<Domain*, Domain*> growthpoint {&staple1_d_1, &staple2_d_2};
        movetype.unassign_domains(staple);
        movetype.set_growthpoint_and_grow_staple(growthpoint, staple);
        delta_e = origami.hybridization_energy(staple1_d_1, staple2_d_2);
        double expected_old_bias {exp(-delta_e) * (5*6)};

        THEN("Old bias matches hand calculated values") {
            REQUIRE(Approx(expected_old_bias) == movetype.m_bias);
        }
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
        CTCBScaffoldRegrowthMCMovetype movetype {origami, random_gens,
                ideal_random_walks, params};

        // Scaffold domains to be regrown
        vector<Domain*> scaffold_domains {&scaffold_d_2, &scaffold_d_3, &scaffold_d_4};

        // The correctness of this is checked in another scneario
        movetype.m_constraintpoints.calculate_constraintpoints(scaffold_domains);

        // Find staples to be regrown (correctness not directly checked)
        set<int> staples {movetype.m_constraintpoints.staples_to_be_regrown()};

        // Unassign scaffold domains to be regrown
        movetype.unassign_domains({&scaffold_d_3, &scaffold_d_4});
        for (auto c_i: staples) {
            movetype.unassign_domains(origami.get_chain(c_i));
        }

        // Set staple 2 domain 2
        vector<Domain*> staple2 {&staple2_d_1, &staple2_d_2};
        double delta_e {origami.set_domain_config(staple2_d_2, {1, 0, 0}, {0, -1, 0})};
        movetype.m_constraintpoints.update_endpoints(&staple2_d_2);
        movetype.m_bias *= exp(-delta_e);

        // Grow staple 2 domain 1
        vector<pair<VectorThree, VectorThree>> configs {};
        vector<double> bfactors {};
        VectorThree p_prev {staple2_d_2.m_pos};
        movetype.calc_biases(staple2_d_1, p_prev, configs, bfactors);
        movetype.calc_bias(bfactors, &staple1_d_1, configs, p_prev, staple2);
        origami.set_domain_config(staple2_d_1, {1, 1, 0}, {0, -1, 0});
        movetype.m_constraintpoints.update_endpoints(&staple2_d_1);

        // Grow staple 3
        configs = {};
        bfactors = {};
        p_prev = scaffold_d_2.m_pos;
        movetype.calc_biases(scaffold_d_3, p_prev, configs, bfactors);
        movetype.calc_bias(bfactors, &scaffold_d_3, configs, p_prev, scaffold_domains);
        origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 1, 0});
        movetype.m_constraintpoints.update_endpoints(&scaffold_d_3);

        // Grow staple 4
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

        THEN("Trial bias matches hand calculated value") {
            REQUIRE(Approx(expected_new_bias) == movetype.m_bias);
        }

        // Regrow old
        movetype.setup_for_regrow_old();
        movetype.m_constraintpoints.reset_active_endpoints();
        movetype.unassign_domains({&scaffold_d_3, &scaffold_d_4});
        for (auto c_i: staples) {
            movetype.unassign_domains(origami.get_chain(c_i));
        }
        movetype.grow_staple_and_update_endpoints(scaffold_domains[0]);
        movetype.grow_chain(scaffold_domains);

        // Old and new bias happen to be the same in this case (maybe do another)
        THEN("Old bias matches hand calculated value") {
            REQUIRE(Approx(expected_new_bias) == movetype.m_bias);
        }
    }
}

SCENARIO("Connector staples are correctly identified") {
    double temp {300};
    double cation_M {1};
    RandomGens random_gens {};
    IdealRandomWalks ideal_random_walks {};
    InputParameters params {};

    // Using a four domain scaffold system for all tests here
    OrigamiSystem origami {setup_four_domain_scaffold_origami(temp, cation_M)};

    origami.add_chain(1);
    origami.add_chain(1);
    origami.add_chain(1);

    // Easy reference
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

    // Also to all domains of a staple
    vector<Domain*> staple11 {&staple11_d_1, &staple11_d_2};
    vector<Domain*> staple12 {&staple12_d_1, &staple12_d_2};
    vector<Domain*> staple13 {&staple13_d_1, &staple13_d_2};

    // Set scaffold configuration
    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, -1, 0});
    origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 1, 0});
    origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 1, 0});
    origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, -1, 0});

    IdentityMCMovetype movetype {origami, random_gens,
            ideal_random_walks, params};

    // These are drawn out in my notebook
    WHEN("Case 1: one staple binding the terminal scaffold domains") {
        origami.set_domain_config(staple11_d_1, {0, 1, 0}, {0, 1, 0});
        origami.set_domain_config(staple11_d_2, {0, 0, 0}, {0, 1, 0});

        THEN("Staple not a connector") {
            REQUIRE(movetype.staple_is_connector(staple11) == false);
        }
    }

    WHEN("Case 2: One staple bound only to another staple") {
        origami.set_domain_config(staple11_d_1, {1, 1, 0}, {0, -1, 0});
        origami.set_domain_config(staple11_d_2, {2, 1, 0}, {0, 1, 0});
        origami.set_domain_config(staple12_d_1, {2, 1, 0}, {0, -1, 0});
        origami.set_domain_config(staple12_d_2, {3, 1, 0}, {0, 1, 0});

        THEN("Staple bound to scaffold and other staple is a connector") {
            REQUIRE(movetype.staple_is_connector(staple11) == true);
        }
    }

    WHEN("Case 3: One staple bound only to another staple bound to another "
            "bound to the scaffold") {
        origami.set_domain_config(staple11_d_1, {1, 1, 0}, {0, -1, 0});
        origami.set_domain_config(staple11_d_2, {2, 1, 0}, {0, 1, 0});
        origami.set_domain_config(staple12_d_1, {2, 1, 0}, {0, -1, 0});
        origami.set_domain_config(staple12_d_2, {3, 1, 0}, {0, 1, 0});
        origami.set_domain_config(staple13_d_1, {3, 1, 0}, {0, -1, 0});
        origami.set_domain_config(staple13_d_2, {4, 1, 0}, {0, 1, 0});

        THEN("Staple between staple bound to scaffold and other staple is a "
                "connector") {
            REQUIRE(movetype.staple_is_connector(staple12) == true);
        }
    }

    WHEN("Case 4: One staple bound to two others which are both bound to scaffold") {
        origami.set_domain_config(staple11_d_1, {1, 1, 0}, {0, -1, 0});
        origami.set_domain_config(staple11_d_2, {2, 1, 0}, {0, 1, 0});
        origami.set_domain_config(staple12_d_1, {2, 1, 0}, {0, -1, 0});
        origami.set_domain_config(staple12_d_2, {2, 0, 0}, {0, 1, 0});
        origami.set_domain_config(staple13_d_1, {2, 0, 0}, {0, -1, 0});
        origami.set_domain_config(staple13_d_2, {1, 0, 0}, {0, -1, 0});

        THEN("Staple in the middle is not a connector") {
            REQUIRE(movetype.staple_is_connector(staple12) == false);
        }
    }
}

void simple_move_iteration(OrigamiSystem& origami, MCMovetype& movetype) {
    // Attempt move and reset if rejected
    Chains original_chains {origami.chains()};
    double original_energy {origami.energy()};
    bool accepted {movetype.attempt_move()};

    // The movetype resets to new configuration if accepted
    if (not accepted) {
        movetype.reset_origami();
        Chains new_chains {origami.chains()};
        double new_energy {origami.energy()};
        THEN("Rejected attempts reset correctly") {
            REQUIRE(compare_vector_contents(original_chains, new_chains));
            REQUIRE(Approx(original_energy) == new_energy);
        }
    }
}

SCENARIO("Moves are fully run and origami system reset each time") {
    double temp {300};
    double cation_M {1};
    RandomGens random_gens {};
    IdealRandomWalks ideal_random_walks {};
    InputParameters params {};

    OrigamiSystem origami {setup_four_domain_scaffold_origami(temp, cation_M)};
    origami.add_chain(1);
    origami.add_chain(2);

    // Easy reference
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

    // Set initial configuration
    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, -1, 0});
    origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 1, 0});
    origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 1, 0});
    origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, -1, 0});
    origami.set_domain_config(staple1_d_1, {0, 1, 0}, {0, 1, 0});
    origami.set_domain_config(staple1_d_2, {0, 0, 0}, {0, 1, 0});
    origami.set_domain_config(staple2_d_1, {1, 0, 0}, {0, -1, 0});
    origami.set_domain_config(staple2_d_2, {1, 1, 0}, {0, -1, 0});

    // Each movetype tried ten times
    WHEN("Met staple exchanges are carried out") {
        for (int i {0}; i != 10; i++) {
            MetStapleExchangeMCMovetype movetype {origami,
                    random_gens, ideal_random_walks, params};
            simple_move_iteration(origami, movetype);
        }
    }

    WHEN("CB staple regrowths are carried out") {
        for (int i {0}; i != 10; i++) {
            CBStapleRegrowthMCMovetype movetype {origami,
                    random_gens, ideal_random_walks, params};
            simple_move_iteration(origami, movetype);
        }
    }

    WHEN("CTCB scaffold regrowths are carried out") {
        for (int i {0}; i != 10; i++) {
            CTCBScaffoldRegrowthMCMovetype movetype {origami, random_gens,
                    ideal_random_walks, params};
            simple_move_iteration(origami, movetype);
        }
    }
}
