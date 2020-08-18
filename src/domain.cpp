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
        VectorThree ndr_3 {cd_j4.m_pos - cd_j3.m_pos};
        bool pair1_stacked {false};
        if (ndr_1 != this->m_ore and ndr_1 != -this->m_ore) {
            pair1_stacked = this->check_twist_constraint(ndr_1, cd_j2);
        }
        if (pair1_stacked) {
            bool pair2_stacked {false};
            if (ndr_3 != cd_j3.m_ore and ndr_3 != -cd_j3.m_ore) {
                pair2_stacked = cd_j3.check_twist_constraint(ndr_3, cd_j4);
            }
            if (pair2_stacked) {
                if (this->m_d > cd_j2.m_d) {
                    ndr_1 = -ndr_1;
                }
                if (not cd_k1.check_twist_constraint(ndr_1, cd_k2)) {
                    kink_constraint_obeyed = false;
                }
            }
        }
    }
    return kink_constraint_obeyed;
}

bool check_pair_stacked(Domain* cd_1, Domain* cd_2) {

    bool stacked {false};
    if (cd_1->m_d > cd_2->m_d) {
        Domain* hold {cd_1};
        cd_1 = cd_2;
        cd_2 = hold;
    }
    VectorThree ndr {cd_2->m_pos - cd_1->m_pos};
    if (ndr != cd_1->m_ore and ndr != -cd_1->m_ore) {
        stacked = cd_1->check_twist_constraint(ndr, *cd_2);
    }

    return stacked;
}
} // namespace domainContainer
