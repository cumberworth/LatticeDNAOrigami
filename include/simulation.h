// simulation.h

#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <memory>
#include <map>

#include "origami_system.h"
#include "movetypes.h"
#include "files.h"

using std::vector;
using std::unique_ptr;
using std::map;

using namespace Origami;
using namespace Movetypes;
using namespace Files;

namespace Simulation {

    class GCMCSimulation {
        public:
            GCMCSimulation(
                    OrigamiSystem origami_system,
                    OrigamiTrajOutputFile output_file,
                    vector<MovetypeConstructor> movetype_constructors,
                    vector<double> movetype_probs);

            void run(int steps, int logging_freq, int center_freq);

        private:

            unique_ptr<MCMovetype> select_movetype();
            void write_log_entry(int step);

            // Big things
            OrigamiSystem m_origami_system;
            OrigamiTrajOutputFile m_output_file;

            // Movetype information
            vector<MovetypeConstructor> m_movetype_constructors;
            vector<double> m_cumulative_probs;

    };

}

#endif // SIMULATION_H

