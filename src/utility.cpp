// utility.cpp

#include <iostream>
#include <vector>
#include <memory>
#include <algorithm>

#include "utility.h"

namespace Utility {

    using std::cout;
    using std::vector;
    using std::unique_ptr;

    int index(vector<int> container, int element) {
        for (size_t i {0}; i != container.size(); i++) {
            if (container[i] == element) {
                return i;
            }
            else continue;
        }
        throw NoElement {};
    }

    VectorThree VectorThree::operator-() {
        VectorThree neg {};
        for (size_t i {0}; i != 3; i++) {
            neg[i] = -m_container[i];
        }
        return neg;
    }

    VectorThree VectorThree::operator+(const VectorThree& v_2) const {
        VectorThree sum;
        for (unsigned int i {0}; i != 3; i++) {
            sum[i] = m_container[i] + v_2.at(i);
        }
        return sum;
    }

    VectorThree VectorThree::operator-(const VectorThree& v_2) const {
        VectorThree diff;
        for (unsigned int i {0}; i != 3; i++) {
            diff[i] = m_container[i] - v_2.at(i);
        }
        return diff;
    }

    bool VectorThree::operator!=(const VectorThree& v_2) const {
        for (unsigned int i {0}; i != 3; i++) {
            if (m_container[i] != v_2.at(i)) {
                return true;
            }
        }
        return false;
    }

    bool operator==(const VectorThree& v1, const VectorThree& v2) {
        for (unsigned int i {0}; i != 3; i++) {
            if (v1.at(i) != v2.at(i)) {
                return false;
            }
        }
        return true;
    }

    VectorThree VectorThree::rotate_half(VectorThree axis) {
        VectorThree rot = *this;
        axis = axis.absolute();
        if (axis == xhat) {
            rot[1] = -rot[1];
            rot[2] = -rot[2];
        }
        else if (axis == yhat) {
            rot[0] = -rot[0];
            rot[2] = -rot[2];
        }
        else if (axis == zhat) {
            rot[0] = -rot[0];
            rot[1] = -rot[1];
        }
//        else {
//           throw SimulationMisuse {};
//      }
        return rot;
    }

    VectorThree VectorThree::rotate(VectorThree origin, VectorThree axis,
            int turns) {

        VectorThree rot = *this;
        if (turns == 0) {
            return rot;
        }
        rot = rot - origin;
        if (turns % 2 == 0) {
            if (axis == xhat) {
                rot[1] = -rot[1];
                rot[2] = -rot[2];
            }
            else if (axis == yhat) {
                rot[0] = -rot[0];
                rot[2] = -rot[2];
            }
            else if (axis == zhat) {
                rot[0] = -rot[0];
                rot[1] = -rot[1];
            }
            else {
                throw SimulationMisuse {};
            }
        }
        else {
            int dir;
            bool odd_turns_even {(turns - 1)/2 % 2 == 0};
            bool turns_neg {turns < 0};
            if (not turns_neg and odd_turns_even) {
                dir = 1;
            }
            else {
                dir = -1;
            }

            if (axis == xhat) {
                std::swap(rot[1], rot[2]);
                rot[1] *= -dir;
                rot[2] *= dir;
            }
            else if (axis == yhat) {
                std::swap(rot[0], rot[2]);
                rot[0] *= -dir;
                rot[2] *= dir;
            }
            else if (axis == zhat) {
                std::swap(rot[0], rot[1]);
                rot[0] *= -dir;
                rot[1] *= dir;
            }
            else {
                throw SimulationMisuse {};
            }
        }

        return rot + origin;
    }

    int VectorThree::sum() {
        return this->at(0) + this->at(1) + this->at(2);
    }

    int VectorThree::abssum() {
        return abs(m_container[0]) + abs(m_container[1]) + abs(m_container[2]);
    }

    VectorThree VectorThree::absolute() {
        return {abs(m_container[0]), abs(m_container[1]), abs(m_container[2])};
    }

    VectorThree VectorThree::sort() {
        array<int, 3> sorted_container {m_container};
        std::sort(sorted_container.begin(), sorted_container.end());
        int x {sorted_container[0]};
        int y {sorted_container[1]};
        int z {sorted_container[2]};
        return {VectorThree {x, y, z}};
    }
}
