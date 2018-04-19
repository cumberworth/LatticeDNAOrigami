// domain.h

#ifndef DOMAIN_H
#define DOMAIN_H

#include "utility.h"

namespace domainContainer{

    using utility::Occupancy;
    using utility::VectorThree;

    class Domain {
        // DNA origami binding domain
        public:

            // Standard methods
            Domain(int c, int c_ident, int d, int d_ident, int c_length): m_c {c},
                    m_c_ident {c_ident}, m_d {d}, m_d_ident {d_ident},
                    m_c_length {c_length} {};
            virtual ~Domain() = default;

            // Immutable attributes
            const int m_c; // Unique index of the associated chain
            const int m_c_ident; // Identity of the chain the associated chain
            const int m_d; // Domain index
            const int m_d_ident; // Domain identity
            const int m_c_length; // Associated chain length
            Domain* m_forward_domain {nullptr}; // Domain in ?' direction
            Domain* m_backward_domain {nullptr}; // Domain in ?' direction

            // State attributes
            VectorThree m_pos {}; // Position vector
            VectorThree m_ore {}; // Orientation vector
            Occupancy m_state {Occupancy::unassigned}; // Binding state
            Domain* m_bound_domain {nullptr}; // Pointer to bound domain

            // Get next domain in chain
            Domain* operator+(int increment);

            // Constraint checkers
            virtual bool check_twist_constraint(VectorThree ndr, Domain& cd_j) = 0;
            virtual bool check_kink_constraint(VectorThree ndr, Domain& cd_j) = 0;
    };

    class SixteenDomain: public Domain {
        using Domain::Domain;
        public:

            // Constraint checkers
            bool check_twist_constraint(VectorThree ndr, Domain& cd_j);
            bool check_kink_constraint(VectorThree ndr, Domain& cd_j);
    };

}

#endif // DOMAIN_H
