// utility.cpp

#include "utility.h"

using namespace Utility;

template<typename Container_T, typename Element_T>
int index(Container_T container, Element_T element) {
    for (int i {0}; i != container.size; i++) {
        if (container[i] == element) {
            return i;
        }
        else continue;
    }
    throw NoElement {};
}

VectorThree VectorThree::operator+(const VectorThree& v_2) {
    VectorThree sum;
    for (int i {0}; i != size(); i++) {
        sum[i] = (*this)[i] + v_2[i];
    }
    return sum;
}

VectorThree VectorThree::operator-(const VectorThree& v_2) {
    VectorThree diff;
    for (int i {0}; i != size(); i++) {
        diff[i] = v_2[i] - (*this)[i];
    }
    return diff;
}

VectorThree VectorThree::operator-() {
    VectorThree neg;
    for (auto i: *this) {
        neg[i] = -i;
    }
    return neg;
}

VectorThree VectorThree::rotate_half(VectorThree axis) {
    VectorThree rot = *this;
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
    return rot;
}
