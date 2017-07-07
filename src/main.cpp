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
#include "order_params.h"
#include "enumerate.h"

int main(int argc, char* argv[]) {

    using std::cout;
    using namespace Parser;
    using namespace Utility;
    using namespace Files;
    using namespace Origami;
    using namespace DomainContainer;
    using namespace Simulation;
    using namespace ConstantTemp;
    using namespace Annealing;
    using namespace PTMC;
    using namespace US;
    using namespace Movetypes;
    using namespace Enumerator;

    InputParameters params {argc, argv};
    OrigamiSystem* origami {setup_origami(params)};

    // Enumerate or simulate
    if (params.m_simulation_type == "enumerate") {
        enumerate_main(*origami, params);
    }
    else {

        // Select simulation type
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
        else {
            cout << "No such simulation type.\n";
            std::exit(1);
        }

        sim->run();
        delete sim;
        delete origami;
    }
}
