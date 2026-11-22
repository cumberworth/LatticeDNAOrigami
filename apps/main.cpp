// main.cpp

#include <iostream>
#include <memory>

#include "LatticeDNAOrigami/annealing_simulation.hpp"
#include "LatticeDNAOrigami/bias_functions.hpp"
#include "LatticeDNAOrigami/constant_temp_simulation.hpp"
#include "LatticeDNAOrigami/domain.hpp"
#include "LatticeDNAOrigami/enumerate.hpp"
#include "LatticeDNAOrigami/files.hpp"
#include "LatticeDNAOrigami/nearest_neighbour.hpp"
#include "LatticeDNAOrigami/order_params.hpp"
#include "LatticeDNAOrigami/parser.hpp"
#include "LatticeDNAOrigami/ptmc_simulation.hpp"
#include "LatticeDNAOrigami/us_simulation.hpp"
#include "LatticeDNAOrigami/version.hpp"

int main(int argc, char* argv[]) {
    try {

        using std::cout;

        parser::InputParameters params {argc, argv};
        // Declared before the simulation, so it is destroyed after it
        std::unique_ptr<origami::OrigamiSystem> origami {
                origami::setup_origami(params)};
        orderParams::SystemOrderParams& ops {
                origami->get_system_order_params()};
        biasFunctions::SystemBiases& biases {origami->get_system_biases()};

        cout << "Git commit hash: " << GIT_COMMIT << "\n";

        // Enumerate or simulate
        if (params.m_simulation_type == "enumerate") {
            cout << "Running enumeration\n";
            enumerator::enumerate_main(*origami, ops, biases, params);
        }
        else {

            // Select simulation type
            std::unique_ptr<simulation::GCMCSimulation> sim;
            if (params.m_simulation_type == "constant_temp") {
                cout << "Running serial constant temperature simulation\n";
                sim = std::make_unique<constantTemp::ConstantTGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else if (params.m_simulation_type == "annealing") {
                cout << "Running serial annealing simulation\n";
                sim = std::make_unique<annealing::AnnealingGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else if (params.m_simulation_type == "t_parallel_tempering") {
                cout << "Running T parallel tempering simulation\n";
                sim = std::make_unique<ptmc::TPTGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else if (params.m_simulation_type == "ut_parallel_tempering") {
                cout << "Running uT parallel tempering simulation\n";
                sim = std::make_unique<ptmc::UTPTGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else if (params.m_simulation_type == "hut_parallel_tempering") {
                cout << "Running HuT parallel tempering simulation\n";
                sim = std::make_unique<ptmc::HUTPTGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else if (params.m_simulation_type == "st_parallel_tempering") {
                cout << "Running ST parallel tempering simulation\n";
                sim = std::make_unique<ptmc::STPTGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else if (params.m_simulation_type == "2d_parallel_tempering") {
                cout << "Running 2D (T and stacking) parallel tempering "
                        "simulation\n";
                sim = std::make_unique<ptmc::TwoDPTGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else if (params.m_simulation_type == "umbrella_sampling") {
                cout << "Running single window US simulation\n";
                sim = std::make_unique<us::SimpleUSGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else if (params.m_simulation_type == "mw_umbrella_sampling") {
                cout << "Running multi-window US simulation\n";
                sim = std::make_unique<us::MWUSGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else if (params.m_simulation_type == "ptmw_umbrella_sampling") {
                cout << "Running parallel tempering multi-window US "
                        "simulation\n";
                sim = std::make_unique<us::PTMWUSGCMCSimulation>(
                        *origami, ops, biases, params);
            }
            else {
                throw utility::NotImplemented {"No such simulation type"};
            }

            sim->run();
        }
    } catch (const std::exception& e) {
        std::cout << std::endl;
        std::cout << "An exception occurred during the run" << std::endl;
        std::cout << std::endl;
        std::cerr << e.what() << std::endl;
        std::cout << std::endl;
        std::cout << "Ending run unsuccesfully" << std::endl;
        std::cout << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
