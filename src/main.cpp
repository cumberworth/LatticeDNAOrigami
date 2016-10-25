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
    InputParameters input_parameters {argc, argv};

    // Setup origami system
    OrigamiInputFile origami_input {input_parameters.m_origami_input_filename};
    vector<vector<int>> identities {origami_input.m_identities};
    vector<vector<string>> sequences {origami_input.m_sequences};
    vector<Chain> configs {origami_input.m_chains};
    OrigamiSystem origami {
            identities,
            sequences,
            configs,
            input_parameters.m_temp,
            input_parameters.m_staple_M,
            input_parameters.m_cation_M,
            input_parameters.m_lattice_site_volume,
            input_parameters.m_cyclic};

    // Setup simulation
    OrigamiTrajOutputFile config_out {input_parameters.m_configs_output_filename,
            input_parameters.m_configs_output_freq, origami};
    OrigamiCountsOutputFile counts_out {input_parameters.m_counts_output_filename,
            input_parameters.m_counts_output_freq, origami};
    vector<OrigamiOutputFile*> outs {&config_out, &counts_out};
    GCMCSimulation sim {origami,
        outs,
        input_parameters.m_movetype_constructors,
        input_parameters.m_movetype_probs};

    if (input_parameters.m_simulation_type == "constant_temp") {
        sim.run_constant_temp(input_parameters.m_steps,
                input_parameters.m_logging_freq,
                input_parameters.m_centering_freq);
    }
    else if (input_parameters.m_simulation_type == "annealing") {
        sim.run_annealing(input_parameters.m_steps_per_temp,
                input_parameters.m_max_temp,
                input_parameters.m_min_temp,
                input_parameters.m_temp_interval,
                input_parameters.m_logging_freq,
                input_parameters.m_centering_freq);
    }
}
