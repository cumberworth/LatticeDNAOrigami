// constant_temp_simulation.h

#ifndef CONSTANT_TEMP_SIMULATION_H
#define CONSTANT_TEMP_SIMULATION_H

#include "simulation.h"

using namespace Simulation;

namespace ConstantTemp {

    class ConstantTGCMCSimulation: public GCMCSimulation {
        public:
            ConstantTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            void run() {simulate(m_steps);}
        private:
            void update_internal(long long int) {};
            long int m_steps;
    };
}

#endif // CONSTANT_TEMP_SIMULATION_H
