#include <catch.hpp>
#include <math.h>

#define private public
#define protected public

#include "parser.h"
#include "files.h"
#include "simulation.h"
#include "test_origami_system.h"

using std::cout;
using std::log;

using namespace Parser;
using namespace Utility;
using namespace Files;
using namespace Origami;
using namespace Simulation;
using namespace TestOrigamiSystem;

// Consider writing some tests for the other simulation methods
//    WHEN("T annealing simulation is run") {
//        //check that t actually increments correctly
//       //check start and end points
//  }


SCENARIO("PTGCMC methods are run", "[!hide][mpi]") {
    double temp {330};
    double cation_M {1};
    InputParameters params {};
    params.m_steps = 1;
    params.m_exchange_interval = 1;
    params.m_num_reps = 2;
    params.m_temps = {330, 340};
    params.m_bias_mults = {1, 1};
    params.m_constant_staple_M = false;
    params.m_staple_M = 1e-6;
    params.m_lattice_site_volume = 4e-28;
    params.m_temp_for_staple_u = 330;
    params.m_chem_pot_mults = {1, 1.5};

    OrigamiSystem origami {setup_two_domain_scaffold_origami(temp, cation_M)};
    origami.add_chain(1);
    SystemBiases system_bias {params, origami};

    PTGCMCSimulation sim {origami, system_bias, params};

    // Easy reference
    Domain& scaffold_d_1 {*origami.get_domain(0, 0)};
    Domain& scaffold_d_2 {*origami.get_domain(0, 1)};
    Domain& staple1_d_1 {*origami.get_domain(1, 0)};
    Domain& staple1_d_2 {*origami.get_domain(1, 1)};


    origami.set_domain_config(scaffold_d_1, {0, 0, 0}, {0, 0, 1});
    origami.set_domain_config(scaffold_d_2, {1, 0, 0}, {0, 0, -1});

    GIVEN("Two replicas with different num staples") {
        if (sim.m_rank == 0) {
            origami.set_domain_config(staple1_d_1, {0, 0, 0}, {0, 0, -1});
            origami.set_domain_config(staple1_d_2, {1, 0, 0}, {0, 0, 1});
        }
        else if (sim.m_rank == 1) {
            origami.delete_chain(1);
        }

        WHEN("PTMC exchange probability is calculated") {
            // Send enthalpy, bias, and num staples from rep 1 to rep 0
            ThermoOfHybrid DH_DS {origami.enthalpy_and_entropy()};
            vector<double> enthalpies {DH_DS.enthalpy};
            vector<double> biases {0};
            vector<int> staples {1};
            if (sim.m_rank == 1) {
                sim.slave_send_and_recieve(0);
            }
            else if (sim.m_rank == 0) {
                sim.master_receive(0, enthalpies, biases, staples);
                sim.master_send(0);
            }

            if (sim.m_rank == 0) {

                // Calculate expected p accept
                double temp1 {330};
                double temp2 {340};
                // mol/L * num/mol * L/lattice sites
                double staple_u1 {temp1 * log(1e-6 * Utility::NA *
                        params.m_lattice_site_volume / 1e-3)};
                double staple_u2 {staple_u1 * 1.5};
                double enthalpy1 {origami.hybridization_enthalpy(scaffold_d_1,
                        staple1_d_1) + origami.hybridization_enthalpy(scaffold_d_2,
                        staple1_d_2)};
                double enthalpy2 {0};
                double DB {1/temp2 - 1/temp1};
                double DH {enthalpy2 - enthalpy1};
                double DBU {staple_u2/temp2 - staple_u1/temp1};
                int DN {-1};
                double exp_p_accept {exp(DB*DH - DBU*DN)};

                double p_accept {sim.calc_acceptance_p(
                        sim.m_temps[0], sim.m_temps[1],
                        sim.m_staple_us[0], sim.m_staple_us[1], 
                        enthalpies[0], enthalpies[1],
                        biases[0], biases[1], staples[0], staples[1])};

                THEN("Value matches expected") {
                   REQUIRE(p_accept == exp_p_accept);
                }
            }
        }
    }
}
