// orientation_movetype.h

#ifndef ORIENTATION_MOVETYPE_H
#define ORIENTATION_MOVETYPE_H

#include "movetypes.h"

namespace movetypes {

    class OrientationRotationMCMovetype:
        public MCMovetype {

        public:
            OrientationRotationMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params);
            OrientationRotationMCMovetype(const
                    OrientationRotationMCMovetype&) = delete;
            OrientationRotationMCMovetype& operator=(const
                    OrientationRotationMCMovetype&) = delete;

            void write_log_summary(ostream* log_stream) override final;

        private:
            bool internal_attempt_move() override;
            void add_external_bias() override final {}
            void add_tracker(bool accepted) override;
    };
}

#endif // ORIENTATION_MOVETYPE_H
