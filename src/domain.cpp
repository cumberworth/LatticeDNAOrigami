// domain.cpp

#include <iostream>

#include "domain.h"

namespace domainContainer {

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

bool HalfTurnDomain::check_twist_constraint(VectorThree ndr, Domain& cd_2) {
    bool twist_constraint_obeyed {true};
    VectorThree ore_1_rotated {m_ore.rotate_half(ndr)};
    if (not(ore_1_rotated == cd_2.m_ore)) {
        twist_constraint_obeyed = false;
    }
    return twist_constraint_obeyed;
}

bool HalfTurnDomain::check_kink_constraint(VectorThree ndr, Domain& cd_2) {
    bool kink_constraint_obeyed {true};
    if (ndr == -m_ore) {
        kink_constraint_obeyed = false;
    }
    else if (ndr == m_ore) {
        if (cd_2.m_ore == -m_ore) {
            kink_constraint_obeyed = false;
        }
    }
    else if (ndr == cd_2.m_ore or ndr == -cd_2.m_ore) {
        kink_constraint_obeyed = false;
    }

    return kink_constraint_obeyed;
}

bool HalfTurnDomain::check_junction_constraint(
        Domain& cd_j2,
        Domain& cd_j3,
        Domain& cd_j4,
        Domain& cd_k1,
        Domain& cd_k2) {

    return true;
}

bool ThreeQuarterTurnDomain::check_twist_constraint(
        VectorThree ndr,
        Domain& cd_2) {
    bool twist_constraint_obeyed {true};
    VectorThree ore_1_rotated {m_ore.rotate(ndr, -1)};
    if (not(ore_1_rotated == cd_2.m_ore)) {
        twist_constraint_obeyed = false;
    }
    return twist_constraint_obeyed;
}

bool ThreeQuarterTurnDomain::check_kink_constraint(
        VectorThree ndr,
        Domain& cd_2) {
    bool kink_constraint_obeyed {true};
    if (ndr == -m_ore) {
        kink_constraint_obeyed = false;
    }
    else if (ndr == m_ore) {
        if (m_ore == cd_2.m_ore or m_ore == -cd_2.m_ore) {
            kink_constraint_obeyed = false;
        }
    }
    else if (ndr == cd_2.m_ore or ndr == -cd_2.m_ore) {
        kink_constraint_obeyed = false;
    }

    return kink_constraint_obeyed;
}

bool ThreeQuarterTurnDomain::check_junction_constraint(
        Domain& cd_j2,
        Domain& cd_j3,
        Domain& cd_j4,
        Domain& cd_k1,
        Domain& cd_k2) {

    bool kink_constraint_obeyed {true};
    VectorThree ndr_k1 {cd_k2.m_pos - cd_k1.m_pos};
    if (ndr_k1 == cd_k1.m_ore) {
        VectorThree ndr_1 {cd_j2.m_pos - this->m_pos};
        if (this->m_d > cd_j2.m_d) {
            ndr_1 = -ndr_1;
        }
        if (not cd_k1.check_twist_constraint(ndr_1, cd_k2)) {
            kink_constraint_obeyed = false;
        }
    }
    return kink_constraint_obeyed;
}
} // namespace domainContainer
