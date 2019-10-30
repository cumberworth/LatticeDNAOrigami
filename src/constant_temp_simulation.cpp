// constant_temp_simulation.cpp

#include <iostream>

#include "constant_temp_simulation.h"

namespace constantTemp {

using std::cout;

ConstantTGCMCSimulation::ConstantTGCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        GCMCSimulation(origami_system, ops, biases, params),
        m_steps {params.m_ct_steps} {

    m_logging_stream = &cout;
    m_output_files = simulation::setup_output_files(
            params,
            params.m_output_filebase,
            m_origami_system,
            m_ops,
            m_biases);
}

vector<double> ConstantTGCMCSimulation::get_energies() { return m_enes; }

vector<int> ConstantTGCMCSimulation::get_staples() { return m_staples; }

vector<int> ConstantTGCMCSimulation::get_domains() { return m_domains; }

void ConstantTGCMCSimulation::update_internal(long long int step) {
    if (m_op_freq != 0 and step % m_op_freq == 0) {
        m_enes.push_back(m_origami_system.energy());
        m_staples.push_back(m_origami_system.num_staples());
        m_domains.push_back(m_origami_system.num_fully_bound_domain_pairs());
    }
}
} // namespace constantTemp
