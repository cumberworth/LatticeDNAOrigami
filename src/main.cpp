// main.cpp

#include <iostream>

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

int main(int argc, char* argv[]) {
    InputParameters params {argc, argv};

    // Setup origami system
    OrigamiInputFile origami_input {params.m_origami_input_filename};
    vector<vector<int>> identities {origami_input.m_identities};
    vector<vector<string>> sequences {origami_input.m_sequences};
    vector<Chain> configs {origami_input.m_chains};
    OrigamiSystem origami {
            identities,
            sequences,
            configs,
            params.m_temp,
            params.m_staple_M,
            params.m_cation_M,
            params.m_lattice_site_volume,
            params.m_cyclic};

    // Setup simulation

    GCMCSimulation* sim;
    if (params.m_simulation_type == "constant_temp") {
        sim = new ConstantTGCMCSimulation {origami, params};
    }
    else if (params.m_simulation_type == "annealing") {
        sim = new AnnealingGCMCSimulation {origami, params};
    }
    else if (params.m_simulation_type == "parallel_tempering") {
        sim = new PTGCMCSimulation {origami, params};
    }
    else {
        cout << "No such simulation type.\n";
        std::exit(1);
    }
    sim->run();
    delete sim;
}
