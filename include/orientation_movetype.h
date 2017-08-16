// orientation_movetype.h

#ifndef ORIENTATION_MOVETYPE_H
#define ORIENTATION_MOVETYPE_H

#include "movetypes.h"

namespace movetypes {

    class OrientationRotationMCMovetype:
        public MCMovetype {

        public:
            using MCMovetype::MCMovetype;

            void write_log_summary(ostream* log_stream) override final;

        private:
            bool internal_attempt_move() override;
            void add_external_bias() override final {}
            void add_tracker(bool accepted) override;
    };
}

#endif // ORIENTATION_MOVETYPE_H
