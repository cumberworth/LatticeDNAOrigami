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

SCENARIO("Distributions are calculated for methods with random variables") {
    double eps {0.01}; // Check values are within 1%
    double temp {300};
    double cation_M {1};
    RandomGens random_gens {};
    IdealRandomWalks ideal_random_walks {};
    InputParameters params {};

    OrigamiSystem origami {setup_two_domain_scaffold_origami(temp, cation_M)};
    SystemBias system_bias {params, origami};
    origami.add_chain(1);

    // Easy reference
    Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
    Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
    Domain& staple_d_1 {*origami.get_domain(1, 0)};
    Domain& staple_d_2 {*origami.get_domain(1, 1)};

    // Just need the shared movetype methods
    IdentityMCMovetype movetype {origami, system_bias, random_gens,
            ideal_random_walks, params};

    GIVEN("Uniform distribution for two domain one staple system domain selection.") {
        unordered_map<pair<int, int>, double> exp_dist {
                {{0, 0}, 0.25},
                {{0, 1}, 0.25},
                {{1, 0}, 0.25},
                {{1, 1}, 0.25}};

        // Estimate distribution from function
        unordered_map<pair<int, int>, double> calc_dist {};
        int num_iters {100000};
        for (int i {0}; i != num_iters; i++) {
            Domain* domain {movetype.select_random_domain()};
            pair<int, int> key {domain->m_d, domain->m_c};
            calc_dist[key] += 1./num_iters;
        }

        THEN("Distributions match (within 1%)") {
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
    }

    GIVEN("Uniform distrubtion of positions around a previous") {
        VectorThree p_prev {3, 4, 1};
        double exp_p {1./6};

        // Estimate distribution from function
        vector<VectorThree> new_pos {
                {4, 4, 1},
                {2, 4, 1},
                {3, 5, 1},
                {3, 3, 1},
                {3, 4, 2},
                {3, 4, 0}};
        unordered_map<VectorThree, double> calc_dist {};
        int num_iters {100000};
        for (int i {0}; i != num_iters; i++) {
            VectorThree pos {movetype.select_random_position(p_prev)};
            calc_dist[pos] += 1./num_iters;
        }

        THEN("Distributions match (within 1%)") {
            bool dists_match = true;
            for (auto pos: new_pos) {
                double calc_p {calc_dist[pos]};
                if (calc_p < (exp_p - eps) or calc_p > (exp_p + eps)) {
                    dists_match = false;
                }
            }
            REQUIRE(dists_match == true);
        }
    }

    GIVEN("A configuration for the origami system and very high T") {
        origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, -1, 0});
        origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, -1, 0});
        origami.set_domain_config(staple_d_1, {0, 0, 0}, {0, 1, 0});
        origami.set_domain_config(staple_d_2, {0, 1, 0}, {0, 1, 0});

        // Over the simulation, all positions of d2 relative to d1 should be
        // equally probable
        // Could add a WHEN here, and in general to the whole test suite
        double exp_p {1./6};
        vector<VectorThree> new_pos {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0},
                {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};

        // Calculate the distribution of accepted positions for scaffold domain 2
        unordered_map<VectorThree, double> calc_dist {};
        int accepted_moves {0};
        int num_iters {100000};
        for (int i {0}; i != num_iters; i++) {
            CTCBScaffoldRegrowthMCMovetype ct_movetype {origami, system_bias, random_gens,
                    ideal_random_walks, params};
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

        THEN("Distributions match (within 1%)") {
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
