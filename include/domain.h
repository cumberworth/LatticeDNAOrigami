// domain.h

#ifndef DOMAIN_H
#define DOMAIN_H

#include "utility.h"

using Utility::VectorThree;
using Utility::Occupancy;

namespace DomainContainer{

    // Container for chain domain indices
    class Domain {
        public:
            int m_c;
            int m_c_ident;
            int m_d;
            int m_d_ident;
            int m_c_length;
            VectorThree m_pos;
            VectorThree m_occ;
            Occupancy m_state;
            Domain* m_bound_domain;

            Domain(int c, int c_ident, int d, int d_ident, int c_length,
                    Domain* forward_domain, Domain* backward_domain): m_c {c},
                    m_c_ident {c_ident}, m_d {d}, m_d_ident {d_ident},
                    m_c_length {c_length}, m_forward_domain {forward_domain},
                    m_backward_domain {backward_domain} {m_state {Occupancy::unassigned;};
            Domain* operator+(int increment);

        private:
            Domain* m_forward_domain;
            Domain* m_backward_domain;
    };

    // Hash needs this, ambiguous to also have method
    inline bool operator==(const Domain& cd_i, const Domain& cd_j) {return (cd_i.c == cd_j.c
            and cd_i.d == cd_j.d);};

}

#endif // DOMAIN_H
