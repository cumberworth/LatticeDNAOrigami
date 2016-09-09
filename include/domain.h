// domain.h

#ifndef DOMAIN_H
#define DOMAIN_H

#include "utility.h"

using namespace Utility;

namespace DomainContainer{

    // Container for chain domain indices
    class Domain {
        public:
            int m_c;
            int m_c_ident;
            int m_d;
            int m_d_ident;
            int m_c_length;
            VectorThree m_pos {};
            VectorThree m_ore {};
            Occupancy m_state {Occupancy::unassigned};
            Domain* m_bound_domain {nullptr};

            Domain(int c, int c_ident, int d, int d_ident, int c_length): m_c {c},
                    m_c_ident {c_ident}, m_d {d}, m_d_ident {d_ident},
                    m_c_length {c_length} {};
            virtual ~Domain() = default;
            Domain* operator+(int increment);

            virtual void check_twist_constraint(VectorThree ndr, Domain& cd_j) = 0;

            Domain* m_forward_domain {nullptr};
            Domain* m_backward_domain {nullptr};
    };

    class SixteenDomain: public Domain {
        using Domain::Domain;
        public:
//            SixteenDomain(int c, int c_ident, int d, int d_ident, int c_length):
//                Domain {c, c_ident, d, d_ident, c_length} {};
            void check_twist_constraint(VectorThree ndr, Domain& cd_j);
    };

}

#endif // DOMAIN_H
