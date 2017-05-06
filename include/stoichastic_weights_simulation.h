// stoichastic_weights_simulation.h

#ifndef STOICHASTIC_WEIGHTS_SIMULATION_H
#define STOICHASTIC_WEIGHTS_SIMULATION_H

#include "simulation.h"

using namespace Simulation;

namespace StoichasticWeights {

    class SWMCSimulation: public GCMCSimulation {
        public:
            SWMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            void run();
        private:
            void update_internal(long long int) {};
            long long int m_steps;
    };
}

#endif // STOICHASTIC_WEIGHTS_SIMULATION_H
