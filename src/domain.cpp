// domain.cpp

#include <iostream>

#include "domain.h"
#include "utility.h"
#include "origami_system.h"

using std::cout;

using namespace DomainContainer;
using namespace Utility;
using namespace Origami;

Domain* Domain::operator+(int incr) {

    // There is probably a better way to do this
    if (incr > 0) {
        incr -= 1;
        if (m_forward_domain == nullptr) {
            throw IndexOutOfRange {};
        }
        Domain* domain {(*m_forward_domain) + incr};
        return domain;
    }
    else if (incr < 0) {
        incr += 1;
        if (m_backward_domain == nullptr) {
            throw IndexOutOfRange {};
        }
        Domain* domain {(*m_backward_domain) + incr};
        return domain;
    }
    else {
        return this;
    }
}

void SixteenDomain::check_twist_constraint(VectorThree ndr, Domain& cd_2) {
    VectorThree ore_1_rotated {m_ore.rotate_half(ndr)};
    if (ore_1_rotated == cd_2.m_ore) {
            return;
    }
    else {
          throw ConstraintViolation {};
    }
}
