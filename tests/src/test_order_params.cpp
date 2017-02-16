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

SCENARIO("Order parameters and bias functions") {
    double temp {330};
    double staple_M {1e-6};
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
        params.m_restraint_pairs.push_back(0);
        params.m_restraint_pairs.push_back(3);
        params.m_min_dist = 1;
        params.m_max_dist = 3;
        params.m_max_bias = 10;
        params.m_bias_mult = 1;
        SystemBias system_bias {params, origami};

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
            REQUIRE(expected_bias == system_bias.calc_bias());
        }
        GIVEN("Scaffold conformation 2") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_3, {1, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_4, {1, 1, 0}, {0, 0, 0});

            double expected_bias {5};
            REQUIRE(expected_bias == system_bias.calc_bias());
        }
        GIVEN("Scaffold conformation 2") {
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_3, {1, 1, 0}, {0, 0, 0});
            origami.set_domain_config(scaffold_d_4, {0, 1, 0}, {0, 0, 0});

            double expected_bias {0};
            REQUIRE(expected_bias == system_bias.calc_bias());
        }
    }
}
