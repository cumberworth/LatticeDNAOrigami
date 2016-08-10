// domain.cpp

#include "domain.h"
#include "utility.h"

using namespace DomainContainer;
using namespace Utility;

Domain* Domain::operator+(int incr) {

    // There is probably a better way to do this
    if (incr > 0) {
        incr -= 1;
        if (m_forward_domain != nullptr) {
            throw IndexOutOfRange {};
        }
        Domain* domain {m_forward_domain + incr};
        return domain;
    }
    else if (incr < 0) {
        incr -= 1;
        if (m_backward_domain != nullptr) {
            throw IndexOutOfRange {};
        }
        Domain* domain {m_backward_domain + incr};
        return domain;
    }
    else {
        return this;
    }
}
