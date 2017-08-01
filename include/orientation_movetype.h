// orientation_movetype.h

#ifndef ORIENTATION_MOVETYPE_H
#define ORIENTATION_MOVETYPE_H

#include "movetypes.h"

namespace movetypes {

    class OrientationRotationMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move(long long int step) override final;

            void write_log_summary(ostream* log_stream) override final;

        private:
            void add_external_bias() override final {}
            void subtract_external_bias() override final {}
    };
}

#endif // ORIENTATION_MOVETYPE_H
