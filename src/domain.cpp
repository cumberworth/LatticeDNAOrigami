// domain.cpp

#include <iostream>

#include "domain.h"

namespace domainContainer {

    using std::cout;

    Domain* Domain::operator+(int incr) {

        // There is probably a better way to do this
        if (incr > 0) {
            incr -= 1;
            if (m_forward_domain == nullptr) {
                return nullptr;
            }
            Domain* domain {(*m_forward_domain) + incr};
            return domain;
        }
        else if (incr < 0) {
            incr += 1;
            if (m_backward_domain == nullptr) {
                return nullptr;
            }
            Domain* domain {(*m_backward_domain) + incr};
            return domain;
        }
        else {
            return this;
        }
    }

    bool SixteenDomain::check_twist_constraint(VectorThree ndr, Domain& cd_2) {
        bool twist_constraint_obeyed {true};
        VectorThree ore_1_rotated {m_ore.rotate_half(ndr)};
        if (not (ore_1_rotated == cd_2.m_ore)) {
            twist_constraint_obeyed = false;
        }
        return twist_constraint_obeyed;
    }
}
