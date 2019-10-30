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
        if (cd_2.m_ore != m_ore) {
            kink_constraint_obeyed = false;
        }
    }
    else if (cd_2.m_ore == ndr or cd_2.m_ore == -ndr) {
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

    bool junction_stacked {true};
    VectorThree ndr_k1 {cd_k2.m_pos - cd_k1.m_pos};

    // Check if kinked pair forms crossover angle
    if (ndr_k1 != cd_k1.m_ore) {
        VectorThree ndr_1 {cd_j2.m_pos - this->m_pos};
        VectorThree ndr_3 {cd_j4.m_pos - cd_j3.m_pos};

        // Check that the kink next domain vector is orthoganal to the
        // junction pair next domain vector and that the junction pair next
        // domain vector are parallel (sign is unimportant here)
        if ((ndr_1 != ndr_k1 and ndr_1 != -ndr_k1) and
            (ndr_1 == ndr_3 or ndr_1 == -ndr_3)) {
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
                    junction_stacked = false;
                }
            }
        }
    }

    return junction_stacked;
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
    //        else if (ndr == m_ore) {
    //            if (cd_2.m_ore != m_ore.rotate(m_ore, 1)) {
    //                kink_constraint_obeyed = false;
    //            }
    //        }
    else if (cd_2.m_ore == ndr or cd_2.m_ore == -ndr) {
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

    bool junction_stacked {true};
    VectorThree ndr_k1 {cd_k2.m_pos - cd_k1.m_pos};

    // Check if kinked pair forms crossover angle
    if (ndr_k1 != cd_k1.m_ore) {
        VectorThree ndr_1 {cd_j2.m_pos - this->m_pos};
        VectorThree ndr_3 {cd_j4.m_pos - cd_j3.m_pos};

        // Check that the kink next domain vector is orthoganal to the
        // junction pair next domain vector and that the junction pair next
        // domain vector are parallel (sign is unimportant here)
        if ((ndr_1 != ndr_k1 and ndr_1 != -ndr_k1) and
            (ndr_1 == ndr_3 or ndr_1 == -ndr_3)) {
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
                    junction_stacked = false;
                }
            }
        }
    }

    return junction_stacked;
}
} // namespace domainContainer
