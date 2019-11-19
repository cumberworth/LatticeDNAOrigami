// annealing_simulation.cpp

#include <iostream>

#include "annealing_simulation.h"

namespace annealing {

using std::cout;

using origami::molarity_to_chempots;

AnnealingGCMCSimulation::AnnealingGCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        GCMCSimulation(origami_system, ops, biases, params),
        m_max_temp {params.m_max_temp},
        m_min_temp {params.m_min_temp},
        m_temp_interval {params.m_temp_interval},
        m_steps_per_temp {params.m_steps_per_temp},
        m_staple_M {params.m_staple_M},
        m_constant_staple_M {params.m_constant_staple_M} {

    if (fmod(m_max_temp - m_min_temp, m_temp_interval) != 0) {
        cout << "Bad temperature interval";
    }
    m_logging_stream = &cout;
    m_output_files = simulation::setup_output_files(
            params,
            params.m_output_filebase,
            m_origami_system,
            m_ops,
            m_biases);
}

void AnnealingGCMCSimulation::run() {
    double temp {m_max_temp};
    long long int step {0};
    while (temp >= m_min_temp) {
        m_origami_system.update_temp(temp);
        if (not m_constant_staple_M) {
            m_origami_system.update_staple_us(temp, 1);
        }
        step += simulate(m_steps_per_temp, step);
        temp -= m_temp_interval;
    }
}
} // namespace annealing
