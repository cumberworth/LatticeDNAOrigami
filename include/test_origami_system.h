// test_origami_system.h

#ifndef TEST_ORIGAMI_SYSTEM_H
#define TEST_ORIGAMI_SYSTEM_H

#include "origami_system.h"

namespace testOrigamiSystem {

    using namespace origami;

    OrigamiSystem setup_two_domain_scaffold_origami(double temp, double cation_M);
    OrigamiSystem setup_four_domain_scaffold_origami(double temp, double cation_M);
    OrigamiSystem setup_snodin_origami(double temp, double cation_M);

}

#endif // TEST_ORIGAMI_SYSTEM_H
