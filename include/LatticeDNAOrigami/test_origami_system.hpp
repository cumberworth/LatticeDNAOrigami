// test_origami_system.h

#ifndef TEST_ORIGAMI_SYSTEM_H
#define TEST_ORIGAMI_SYSTEM_H

#include "LatticeDNAOrigami/origami_system.hpp"

namespace testOrigamiSystem {

using namespace origami;

OrigamiSystem setup_two_domain_scaffold_origami(double temp, double cation_M);
OrigamiSystem setup_four_domain_scaffold_origami(double temp, double cation_M);
OrigamiSystem setup_snodin_origami(double temp, double cation_M);

} // namespace testOrigamiSystem

#endif // TEST_ORIGAMI_SYSTEM_H
