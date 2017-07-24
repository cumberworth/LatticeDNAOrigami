// annealing_simulation.h

#ifndef ANNEALING_SIMULATION_H
#define ANNEALING_SIMULATION_H

#include "simulation.h"

namespace Annealing {

    using namespace Simulation;

    class AnnealingGCMCSimulation: public GCMCSimulation {
        public:
            AnnealingGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            void run();
        private:
            void update_internal(long long int) {};
            double m_max_temp;
            double m_min_temp;
            double m_temp_interval;
            long long int m_steps_per_temp;
    };
}

#endif // ANNEALING_SIMULATION_H
