// constant_temp_simulation.cpp

#include "stoichastic_weights_simulation.h"
#include "cd_scaffold_regrowth.h"

using namespace StoichasticWeights;
using namespace CDScaffoldRegrowth;

SWMCSimulation::SWMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        GCMCSimulation(origami_system, params),
        m_steps {params.m_steps} {

    m_logging_stream = &cout;
    m_output_files = setup_output_files(params, params.m_output_filebase,
            m_origami_system);

    // Setup params for staple sims
    m_params.m_movetype_probs = {0.5, 0.5};
    m_params.m_movetypes = {2, 4};
    //m_params.m_steps = 2000000;
    m_params.m_steps = 10;
    m_params.m_output_filebase = "staples";
    m_params.m_logging_freq = 0;
    m_params.m_configs_output_freq = 0;
    m_params.m_counts_output_freq = 0;
    m_params.m_order_params_output_freq = 0;
    m_params.m_energies_output_freq = 0;
    m_params.m_centering_freq = 0;
}

void SWMCSimulation::run() {
  
    CDScaffoldRegrowthMCMovetype movetype {m_origami_system,
            m_random_gens, m_ideal_random_walks, m_params};
    movetype.attempt_move();
    movetype.reset_origami();
    double enes_mean {movetype.m_enes_mean};
    double staples_mean {movetype.m_staples_mean};
    double domains_mean {movetype.m_domains_mean};
    cout << staples_mean << " " << domains_mean << "\n";
    for (long long int step {0}; step != (m_steps + 1); step ++) {
  
        // Pick movetype and apply
        CDScaffoldRegrowthMCMovetype movetype {m_origami_system,
                m_random_gens, m_ideal_random_walks, m_params};
        // Make this proper, also only print the previously accepted values
        // Also check how long 4 domain system needs to equilibrate the staple values
        // Should be much shorter than 2 million
        bool accepted;
        movetype.m_prev_enes_mean = enes_mean;
//        movetype.m_burnin = 100;
        movetype.m_burnin = 0;
        accepted = movetype.attempt_move();
        if (not accepted) {
            movetype.reset_origami();
        }
        else {
            enes_mean = movetype.m_enes_mean;
            staples_mean = movetype.m_staples_mean;
            domains_mean = movetype.m_domains_mean;
        }
  
        // Center and check constraints
        //if (m_centering_freq != 0 and step % m_centering_freq == 0) {
            m_origami_system.centre();
            m_origami_system.check_all_constraints();
        //}
  
        // Write log entry to standard out
        *m_logging_stream << "Step: " << step << " ";
        *m_logging_stream << "Movetype: " << movetype.m_label() << " ";
        *m_logging_stream << "Staples: " << staples_mean << " ";
        *m_logging_stream << "Domains: " << domains_mean << " ";
        *m_logging_stream << "Accepted: " << accepted << " ";
        *m_logging_stream << "Temp: " << m_origami_system.m_temp << " ";
        *m_logging_stream << "Energy: " << enes_mean << " ";
        *m_logging_stream << "Bias: " << m_origami_system.bias() << " ";
        *m_logging_stream << "\n";

  
        // Update internal simulation variables
        update_internal(step);
  
        // Write system properties to file
        for (auto output_file: m_output_files) {
            if (output_file->m_write_freq != 0 and step % output_file->m_write_freq == 0) {
                output_file->write(step);
            }
        }
    }
}
