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

OrigamiSystem TestOrigamiSystem::setup_two_domain_scaffold_origami(
        double temp,
        double cation_M) {
    // Setups up system with no staples and unassigned scaffold
    double staple_M {1e-3};
    double lattice_site_volume {4e-28};
    bool cyclic {false};
    string system_filename {"data/two_domain.json"};

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

    origami.unassign_domain(*origami.get_domain(0, 0));
    origami.unassign_domain(*origami.get_domain(0, 1));

    return origami;
}

OrigamiSystem TestOrigamiSystem::setup_four_domain_scaffold_origami(
        double temp,
        double cation_M) {
    // Setups up system with no staples and unassigned scaffold
    // Scaffold: 1 2 3 4, staple 1: 1 4, staple 2: 3 2
    double staple_M {1e-3};
    double lattice_site_volume {4e-28};
    bool cyclic {false};
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

    origami.unassign_domain(*origami.get_domain(0, 0));
    origami.unassign_domain(*origami.get_domain(0, 1));
    origami.unassign_domain(*origami.get_domain(0, 2));
    origami.unassign_domain(*origami.get_domain(0, 3));

    return origami;
}

OrigamiSystem TestOrigamiSystem::setup_snodin_origami(
        double temp,
        double cation_M) {

    double staple_M {1e-3};
    double lattice_site_volume {4e-28};
    bool cyclic {false};
    string system_filename {"data/snodin_unbound.json"};

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

    return origami;
}

SCENARIO("Energies of an origami system are calculated") {

    GIVEN("An origami system with 2 scaffold domains and one 2 domain staple") {
        double temp {300};
        double cation_M {1};
        OrigamiSystem origami {setup_two_domain_scaffold_origami(temp, cation_M)};
        origami.add_chain(1);

        // Easy reference
        Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
        Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
        Domain& staple_d_1 {*origami.get_domain(1, 0)};
        Domain& staple_d_2 {*origami.get_domain(1, 1)};

        // Sequences are both 5' to 3', but staple domains ordered 3' to 5'
        string scaffold_seq_1 {"AATTGGCAGTTTCACA"};
        string scaffold_seq_2 {"CTGCCATACACTAATA"};
        string staple_seq_1 {"TGTGAAACTGCCAATT"};
        string staple_seq_2 {"TATTAGTGTATGGCAG"};

        WHEN("Scaffold 1 and staple 1 bind") {

            // Direct energy calculation with NN
            double scaffold_staple_1_ene {calc_unitless_hybridization_energy(
                    scaffold_seq_1, temp, cation_M)};

            // Energy from origami object
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple_d_1, {0, 0, 0},
                    {-1, 0, 0})};

            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == scaffold_staple_1_ene);
            }
        }

        WHEN("Scaffold 2 and staple 2 bind") {
            // Direct energy calculation with NN
            double scaffold_staple_2_ene {calc_unitless_hybridization_energy(
                    scaffold_seq_2, temp, cation_M)};

            // Energy from origami object
            origami.set_domain_config(scaffold_d_2, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple_d_2, {0, 0, 0},
                    {-1, 0, 0})};

            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == scaffold_staple_2_ene);
            }
        }

        WHEN("Scaffold 1 and scaffold 2 bind") {
            // Direct energy calculation with NN
            vector<string> scaffold_1_2_seqs {find_longest_contig_complement(
                    scaffold_seq_1, scaffold_seq_2)};
            double ene_sum {0};
            for (auto seq: scaffold_1_2_seqs) {
                ene_sum += calc_unitless_hybridization_energy(seq, temp, cation_M);
            }
            double scaffold_1_2_ene {ene_sum / scaffold_1_2_seqs.size()};

            // Energy from origami object
            origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(scaffold_d_2, {0, 0, 0},
                    {-1, 0, 0})};

            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == scaffold_1_2_ene);
            }
        }

        WHEN("Staple 1 and staple 2 bind") {
            // Direct energy calculation with NN
            vector<string> staple_1_2_seqs {find_longest_contig_complement(
                    staple_seq_1, staple_seq_2)};
            double ene_sum {0};
            for (auto seq: staple_1_2_seqs) {
                ene_sum += calc_unitless_hybridization_energy(seq, temp, cation_M);
            }
            double staple_1_2_ene {ene_sum / staple_1_2_seqs.size()};

            // Energy from origami object
            origami.set_domain_config(staple_d_1, {0, 0, 0}, {1, 0, 0});
            double delta_e {origami.set_domain_config(staple_d_2, {0, 0, 0},
                    {-1, 0, 0})};

            THEN("The change in energy is consistent") {
                REQUIRE(Approx(delta_e) == staple_1_2_ene);
            }
        }

        WHEN("Staple 1 and staple 1 bind") {
            // Direct energy calculation with NN
            vector<string> staple_1_1_seqs {find_longest_contig_complement(
                    staple_seq_1, staple_seq_1)};
            double ene_sum {0};
            for (auto seq: staple_1_1_seqs) {
                ene_sum += calc_unitless_hybridization_energy(seq, temp, cation_M);
            }
            double staple_1_1_ene {ene_sum / staple_1_1_seqs.size()};

            // Energy from origami object
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
            // Direct energy calculation with NN
            vector<string> staple_2_2_seqs {find_longest_contig_complement(
                    staple_seq_2, staple_seq_2)};
            double ene_sum {0};
            for (auto seq: staple_2_2_seqs) {
                ene_sum += calc_unitless_hybridization_energy(seq, temp, cation_M);
            }
            double staple_2_2_ene {ene_sum / staple_2_2_seqs.size()};

            // Energy from origami object
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

SCENARIO("Helical constraints are checked") {
    // This is not intended to be comprehensive; the enumeration on the two
    // domain system already ensure that the total number of configuerations is equal
    // to that calculated by hand. This will mainly be to check the junction
    // constraints, but some of the two domain configurations will be checked anyways
    double temp {300};
    double cation_M {1};

    // The following system is for checking the ? constraints
    GIVEN("Origami system with 2 16 nt scaffold domans and 1 2 domain staple") {
        OrigamiSystem origami {setup_two_domain_scaffold_origami(temp, cation_M)};
        origami.add_chain(1);

        // Easy reference
        Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
        Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
        Domain& staple_d_1 {*origami.get_domain(1, 0)};
        Domain& staple_d_2 {*origami.get_domain(1, 1)};

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

    // The following system is for checking the ? constraints
    GIVEN("Origami system with 4 16 nt scaffold domains and 2 2 domain staples") {
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

        // Consider making more clear what type of constraints these three when
        // cases cover, and any further cases they don't
        WHEN("The second staple is bound to its domains in the same helix") {
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

            // Why is there a when same helix, and no other cases?
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

        /* This is the test of the four body junction constraint. I was
            concerned that the order in which the staples are added could affect
            whether the conformation was in violation of the constraints (the code
            is rather convoluted to make sure these different cases are covered)
            so I decided to test every possible order the staples could be added
            both correctly and incorrectly.
        */
        WHEN("Two copies of the second staple are bound to their respective domains") {
            origami.add_chain(1);

            // Easy reference
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

                // Why is there an alternative correctly here but not with the others?
                WHEN("Alt correctly") {
                    origami.set_domain_config(scaffold_d_2, {0, 1, 0}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_3, {0, 1, 1}, {0, 0, 1});
                    origami.set_domain_config(scaffold_d_4, {1, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_1, {0, 1, 1}, {0, 0, -1});
                    origami.set_domain_config(staple2_d_2, {0, 1, 0}, {0, 0, -1});
                    origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, 1});
                    origami.set_domain_config(staple12_d_2, {1, 1, 1}, {0, 0, -1});
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

    // The following system is for checking the ? constraints
    GIVEN("Snodin origami system") {

        OrigamiSystem origami {setup_snodin_origami(temp, cation_M)};
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
