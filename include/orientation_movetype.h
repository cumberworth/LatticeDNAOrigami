// orientation_movetype.h

#ifndef ORIENTATION_MOVETYPE_H
#define ORIENTATION_MOVETYPE_H

#include "movetypes.h"

using namespace Movetypes;

namespace OrientationMovetype {

    class OrientationRotationMCMovetype: public MCMovetype {
        public:
            using MCMovetype::MCMovetype;
            bool attempt_move();

            string m_label() {return "OrientationRotationMCMovetype";};

            //HACK
            void update_bias(int sign) {}
    };
}

#endif // ORIENTATION_MOVETYPE_H