// main.cpp

#include <iostream>

#include "parser.h"
#include "nearest_neighbour.h"
#include "domain.h"
#include "files.h"
#include "simulation.h"
#include "order_params.h"

using std::cout;

using namespace Parser;
using namespace Utility;
using namespace Files;
using namespace Origami;
using namespace DomainContainer;
using namespace Simulation;
using namespace Movetypes;
using namespace OrderParams;

int main(int argc, char* argv[]) {
    InputParameters params {argc, argv};

    // Setup origami system
    OrigamiInputFile origami_input {params.m_origami_input_filename};
    vector<vector<int>> identities {origami_input.m_identities};
    vector<vector<string>> sequences {origami_input.m_sequences};
    vector<Chain> configs;
    if (params.m_restart_traj_file == "") {
        configs = origami_input.m_chains;
    }
    else {
        //OrigamiInputFile origami_input {params.m_origami_input_filename};
    }
    bool cyclic {origami_input.m_cyclic};
    double staple_u {molarity_to_chempot(params.m_staple_M,
            params.m_temp_for_staple_u, params.m_lattice_site_volume)};
    staple_u *= params.m_staple_u_mult;
    double volume {chempot_to_volume(staple_u, params.m_temp)};
    OrigamiSystem origami {
            identities,
            sequences,
            configs,
            params.m_temp,
            volume,
            params.m_cation_M,
            staple_u,
            cyclic,
            params.m_energy_filebase};
    SystemBias system_bias {params, origami};

    // Setup simulation

    GCMCSimulation* sim;
    if (params.m_simulation_type == "constant_temp") {
        sim = new ConstantTGCMCSimulation {origami, system_bias, params};
    }
    else if (params.m_simulation_type == "annealing") {
        sim = new AnnealingGCMCSimulation {origami, system_bias, params};
    }
    else if (params.m_simulation_type == "parallel_tempering") {
        sim = new PTGCMCSimulation {origami, system_bias,params};
    }
    else {
        cout << "No such simulation type.\n";
        std::exit(1);
    }
    sim->run();
    delete sim;
}
