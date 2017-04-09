// constant_temp_simulation.cpp

#include "constant_temp_simulation.h"

using namespace ConstantTemp;

ConstantTGCMCSimulation::ConstantTGCMCSimulation(OrigamiSystem& origami_system,
        InputParameters& params) :
        GCMCSimulation(origami_system, params),
        m_steps {params.m_steps} {

    m_logging_stream = &cout;
    m_output_files = setup_output_files(params, params.m_output_filebase,
            m_origami_system);
}
