#include <catch.hpp>
#include <math.h>

#define private public
#define protected public

#include "parser.h"
#include "files.h"
#include "simulation.h"
#include "ptmc_simulation.h"
#include "test_origami_system.h"

using std::cout;
using std::log;

using namespace Parser;
using namespace Utility;
using namespace Files;
using namespace Origami;
using namespace Simulation;
using namespace PTMC;
using namespace TestOrigamiSystem;

// Consider writing some tests for the other simulation methods
//    WHEN("T annealing simulation is run") {
//        //check that t actually increments correctly
//       //check start and end points
//  }

mpi::environment m_env;
mpi::communicator m_world;

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

    UTPTGCMCSimulation sim {origami, params};

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
            vector<long double> enthalpies {DH_DS.enthalpy};
            vector<double> biases {0};
            vector<int> staples {1};
            vector<vector<double>> dependent_qs {{}, {}, {}};
            if (sim.m_rank == 1) {
                sim.slave_send(0);
                sim.slave_receive(0);
            }
            else if (sim.m_rank == 0) {
                sim.master_receive(0, dependent_qs);
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
                double DH {enthalpy2*temp2 - enthalpy1*temp1};
                double DBU {staple_u2/temp2 - staple_u1/temp1};
                int DN {-1};
                double exp_p_accept {exp(DB*DH - DBU*DN)};

                // Build control quantitiy pairs vector
                vector<pair<double, double>> control_q_pairs {};
                for (auto control_q: sim.m_control_qs) { 
                    double q_1 {control_q[0]};
                    double q_2 {control_q[1]};
                    control_q_pairs.push_back({q_1, q_2});
                 }

                // Build dependent quantity pairs vector
                int repi1 {sim.m_q_to_repi[0]};
                int repi2 {sim.m_q_to_repi[1]};
                vector<pair<double, double>> dependent_q_pairs {};
                for (auto dependent_q: dependent_qs) {
                    double q_1 {dependent_q[repi1]};
                    double q_2 {dependent_q[repi2]};
                    dependent_q_pairs.push_back({q_1, q_2});
                }


                double p_accept {sim.calc_acceptance_p(control_q_pairs,
                        dependent_q_pairs)};

                THEN("Value matches expected") {
                   REQUIRE(p_accept == exp_p_accept);
                }
           }
        }

        WHEN("Exchanges made") {
            double temp_1 {330};
            double temp_2 {340};
            double staple_u_1 {molarity_to_chempot(1e-6, 330., 4e-28)};
            double staple_u_2 {staple_u_1 * 1.5};
            if (sim.m_rank == 0) {

                // Perform swap
                sim.m_q_to_repi[0] = 1;
                sim.m_q_to_repi[1] = 0;

                // Send to replica
                sim.master_send(0);
                THEN("Control qs swapped") {
                    REQUIRE(sim.m_replica_control_qs[sim.m_temp_i] == temp_2);
                    REQUIRE(sim.m_replica_control_qs[sim.m_staple_u_i] == staple_u_2);
                }
            }
            else {
                sim.slave_receive(0);
                THEN("Control qs swapped") {
                    REQUIRE(sim.m_replica_control_qs[sim.m_temp_i] == temp_1);
                    REQUIRE(sim.m_replica_control_qs[sim.m_staple_u_i] == staple_u_1);
                }
            }
        }

    }
}
