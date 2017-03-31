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

/* This really not organized very well. This focuses on testing whether the
    biases are as expected via the origami with bias class interface. It does
    not comprehensivly test all order parameters or bias functions.
*/
SCENARIO("Order parameters and bias functions") {
    double temp {330};
    double staple_M {1e-6};
    double cation_M {1};
    double lattice_site_volume {4e-28};
    bool cyclic {false};
    RandomGens random_gens {};
    IdealRandomWalks ideal_random_walks {};
    InputParameters params {};
    string system_filename {"data/four2_domain_loop.json"};
    OrigamiInputFile origami_input {system_filename};
    vector<vector<int>> identities {origami_input.m_identities};
    vector<vector<string>> sequences {origami_input.m_sequences};
    vector<Chain> configs {origami_input.m_chains};
    double staple_u {molarity_to_chempot(staple_M, temp,
            lattice_site_volume)};
    double volume {chempot_to_volume(staple_u, temp)};

    GIVEN("Four domain loop system with dist bias") {
        // Scaffold: 1 2 3 4, staple 1: 1 4, staple 2: 3 2

        // Order param setup
        params.m_distance_bias = true;
        params.m_distance_pairs.push_back(0);
        params.m_distance_pairs.push_back(3);
        params.m_min_dist = 1;
        params.m_max_dist = 3;
        params.m_max_bias = 10;
        params.m_bias_mult = 1;

        OrigamiSystemWithBias origami {
                identities,
                sequences,
                configs,
                temp,
                volume,
                cation_M,
                staple_u,
                cyclic,
                params};

        Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
        Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
        Domain& scaffold_d_3 {*origami.get_domain(0, 2)};
        Domain& scaffold_d_4 {*origami.get_domain(0, 3)};

        origami.unassign_domain(scaffold_d_1);
        origami.unassign_domain(scaffold_d_2);
        origami.unassign_domain(scaffold_d_3);
        origami.unassign_domain(scaffold_d_4);

        GIVEN("Scaffold conformation 1") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_3, {2, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_4, {3, 0, 0}, {0, 0, 0});

            double expected_bias {10};
            REQUIRE(expected_bias == origami.bias());

            origami.unassign_domain(scaffold_d_4);
            expected_bias = 0;
            REQUIRE(expected_bias == origami.bias());
            origami.unassign_domain(scaffold_d_3);
            origami.unassign_domain(scaffold_d_2);
            origami.unassign_domain(scaffold_d_1);
            REQUIRE(expected_bias == origami.bias());
        }
        
        GIVEN("Scaffold conformation 2") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_3, {1, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_4, {1, 1, 0}, {0, 0, 0});

            double expected_bias {5};
            REQUIRE(expected_bias == origami.bias());
        }

        GIVEN("Scaffold conformation 2") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, 0, 0});

            double expected_bias {0};
            REQUIRE(expected_bias == origami.bias());
        }
    }

    GIVEN("Four domain loop system with dist bias") {
        // Scaffold: 1 2 3 4, staple 1: 1 4, staple 2: 3 2
        OrigamiSystemWithBias origami {
                identities,
                sequences,
                configs,
                temp,
                volume,
                cation_M,
                staple_u,
                cyclic,
                params};

        Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
        Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
        Domain& scaffold_d_3 {*origami.get_domain(0, 2)};
        Domain& scaffold_d_4 {*origami.get_domain(0, 3)};

        origami.unassign_domain(scaffold_d_1);
        origami.unassign_domain(scaffold_d_2);
        origami.unassign_domain(scaffold_d_3);
        origami.unassign_domain(scaffold_d_4);

        // Setup bias grid
        unordered_map<vector<int>, double> bias_grid {};
        bias_grid[{0, 0}] = -10;
        bias_grid[{1, 1}] = -5;
        GridBiasFunction* grid_bias_f {
                origami.get_system_biases()->get_grid_bias()};
        SystemOrderParams* system_order_params {
                origami.get_system_order_params()};
        OrderParam* num_bound_domains {&system_order_params->get_num_bound_domains()};
        OrderParam* num_staples {&system_order_params->get_num_staples()};
        vector<OrderParam*> grid_params {num_bound_domains, num_staples};
        grid_bias_f->set_order_params(grid_params);
        grid_bias_f->replace_biases(bias_grid);

        GIVEN("Unbound scaffold") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_3, {2, 0, 0}, {0, 0, 0});

            double expected_delta_e {-10};
            double delta_e {origami.check_domain_constraints(scaffold_d_4, {3, 0, 0},
                    {0, 0, 0})};
            REQUIRE(expected_delta_e == delta_e);

            origami.set_domain_config(scaffold_d_4, {3, 0, 0}, {0, 0, 0});

            double expected_bias {-10};
            REQUIRE(expected_bias == origami.bias());
        }

        GIVEN("Self-bound scaffold") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_3, {2, 0, 0}, {0, 0, 0});
            double hyb_e {origami.hybridization_energy(scaffold_d_3, scaffold_d_4)};

            double expected_bias {-10};
            double delta_e {origami.check_domain_constraints(scaffold_d_4, {1, 0, 0},
                    {0, 0, 0})};
            REQUIRE(expected_bias == delta_e - hyb_e);

            origami.set_domain_config(scaffold_d_4, {1, 0, 0}, {0, 0, 0});

            REQUIRE(expected_bias == origami.bias());
        }

        GIVEN("Staple bound correctly at one domain and the unbound") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, -1});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 0, 1});
            origami.set_domain_config(scaffold_d_3, {2, 0, 0}, {0, 0, -1});
            origami.set_domain_config(scaffold_d_4, {3, 0, 0}, {0, 0, 1});
            origami.add_chain(1);

            Domain& staple_d_1 {*origami.get_domain(1, 0)};
            Domain& staple_d_2 {*origami.get_domain(1, 1)};

            double delta_e {0};
            delta_e = origami.set_domain_config(staple_d_1, {0, 0, 0}, {0, 0, 1});
            double hyb_e {origami.hybridization_energy(scaffold_d_1, staple_d_1)};
            double expected_delta_e {10 + hyb_e};

            REQUIRE(expected_delta_e == delta_e);

            delta_e = origami.set_domain_config(staple_d_2, {0, 2, 0}, {0, 0, 1});
            expected_delta_e = -5;

            REQUIRE(expected_delta_e == delta_e);

            origami.unassign_domain(staple_d_2);
            delta_e = origami.check_domain_constraints(staple_d_2, {0, 2, 0}, {0, 0, 1});
            REQUIRE(expected_delta_e == delta_e);

            origami.set_checked_domain_config(staple_d_2, {0, 2, 0}, {0, 0, 1});
            double expected_bias {-5};
            REQUIRE(expected_bias == origami.bias());

            origami.unassign_domain(staple_d_2);
            origami.unassign_domain(staple_d_1);
            delta_e = origami.check_delete_chain(1);
            expected_delta_e = -10;
            REQUIRE(expected_delta_e == delta_e);
            origami.delete_chain(1);
            expected_bias = expected_delta_e;
            REQUIRE(expected_bias == origami.bias());
        }
    }
}
