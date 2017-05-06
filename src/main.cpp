// main.cpp

#include <iostream>

#include "parser.h"
#include "nearest_neighbour.h"
#include "domain.h"
#include "files.h"
#include "constant_temp_simulation.h"
#include "annealing_simulation.h"
#include "ptmc_simulation.h"
#include "us_simulation.h"
#include "stoichastic_weights_simulation.h"
#include "order_params.h"
#include "enumerate.h"

using std::cout;

using namespace Parser;
using namespace Utility;
using namespace Files;
using namespace Origami;
using namespace DomainContainer;
using namespace Simulation;
using namespace ConstantTemp;
using namespace Annealing;
using namespace StoichasticWeights;
using namespace PTMC;
using namespace US;
using namespace Movetypes;
using namespace Enumerator;

// This is getting a bit long
int main(int argc, char* argv[]) {
    InputParameters params {argc, argv};

    // Setup origami system
    OrigamiInputFile origami_input {params.m_origami_input_filename};
    vector<vector<int>> identities {origami_input.get_identities()};
    vector<vector<string>> sequences {origami_input.get_sequences()};
    vector<Chain> configs = origami_input.get_config();
    if (params.m_restart_traj_file != "") {
        OrigamiTrajInputFile traj_file {params.m_restart_traj_file};
        configs = traj_file.read_config(params.m_restart_step);
    }

    // Calculate chemical potential from specified staple concentration
    double staple_u {molarity_to_chempot(params.m_staple_M,
            params.m_temp_for_staple_u, params.m_lattice_site_volume)};
    staple_u *= params.m_staple_u_mult;
    double volume {chempot_to_volume(staple_u, params.m_temp)};

    OrigamiSystem* origami;
    if (params.m_biases_present) {
        origami = new OrigamiSystemWithBias {identities, sequences, configs,
                volume, staple_u, params};
    }
    else {
        origami = new OrigamiSystem {identities, sequences, configs,
            volume, staple_u, params};
    }

    // Run enumeration
    if (params.m_simulation_type == "enumerate") {
        enumerate_main(*origami, params);
    }
    else {
        // Run simulation
        GCMCSimulation* sim;
        if (params.m_simulation_type == "constant_temp") {
            sim = new ConstantTGCMCSimulation {*origami, params};
        }
        else if (params.m_simulation_type == "annealing") {
            sim = new AnnealingGCMCSimulation {*origami, params};
        }
        else if (params.m_simulation_type == "t_parallel_tempering") {
            sim = new TPTGCMCSimulation {*origami, params};
        }
        else if (params.m_simulation_type == "ut_parallel_tempering") {
            sim = new UTPTGCMCSimulation {*origami, params};
        }
        else if (params.m_simulation_type == "hut_parallel_tempering") {
            sim = new HUTPTGCMCSimulation {*origami, params};
        }
        else if (params.m_simulation_type == "umbrella_sampling") {
            sim = new SimpleUSGCMCSimulation {*origami, params};
        }
        else if (params.m_simulation_type == "mw_umbrella_sampling") {
            sim = new MWUSGCMCSimulation {*origami, params};
        }
        else if (params.m_simulation_type == "stoichastic_weights") {
            sim = new SWMCSimulation {*origami, params};
        }
        else {
            cout << "No such simulation type.\n";
            std::exit(1);
        }
        sim->run();
        delete sim;
        delete origami;
    }
}
