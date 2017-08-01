// constant_temp_simulation.h

#ifndef CONSTANT_TEMP_SIMULATION_H
#define CONSTANT_TEMP_SIMULATION_H

#include "bias_functions.h"
#include "order_params.h"
#include "origami_system.h"
#include "parser.h"
#include "simulation.h"

namespace constantTemp {

    using biasFunctions::SystemBiases;
    using orderParams::SystemOrderParams;
    using origami::OrigamiSystem;
    using parser::InputParameters;
    using simulation::GCMCSimulation;

    class ConstantTGCMCSimulation: public GCMCSimulation {
        public:
            ConstantTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params);
            void run() {simulate(m_steps);}
            vector<double> get_energies();
            vector<int> get_staples();
            vector<int> get_domains();

            int m_op_freq {0};
        private:
            void update_internal(long long int step);
            long long int m_steps;
            vector<double> m_enes {};
            vector<int> m_staples {};
            vector<int> m_domains {};
    };
}

#endif // CONSTANT_TEMP_SIMULATION_H
