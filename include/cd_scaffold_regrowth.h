// cd_scaffold_regrowth.h

#ifndef CD_SCAFFOLD_REGROWTH_H
#define CD_SCAFFOLD_REGROWTH_H

#include "movetypes.h"
#include "constant_temp_simulation.h"

using namespace Movetypes;
using namespace ConstantTemp;

namespace CDScaffoldRegrowth {

    class CDScaffoldRegrowthMCMovetype: public MCMovetype {
        public:
            CDScaffoldRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    InputParameters& params);
            ~CDScaffoldRegrowthMCMovetype();
            bool attempt_move();

            string m_label() {return "CDScaffoldRegrowthMCMovetype";}

            double m_prev_enes_mean {0};
            double m_enes_mean;
            double m_staples_mean;
            double m_domains_mean;
            int m_burnin {0};

        private:
            double m_delta_e;
            ConstantTGCMCSimulation* m_staple_sim;
    };
}

#endif // CD_SCAFFOLD_REGROWTH_H
