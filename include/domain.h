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
            virtual bool check_junction_constraint(
                    Domain& cd_j2,
                    Domain& cd_j3,
                    Domain& cd_j4,
                    Domain& cd_k1,
                    Domain& cd_k2
                    ) = 0;
    };

    class HalfTurnDomain: public Domain {
        using Domain::Domain;
        public:

            // Constraint checkers
            bool check_twist_constraint(VectorThree ndr, Domain& cd_j);
            bool check_kink_constraint(VectorThree ndr, Domain& cd_j);

            /** Check four body junction stacking rules
              *
              * If the next domain vectors of the first and third pairs are
              * parallel to each other, antiparallel to the kink pair next
              * domain vector, and the kink pair is not agreeing with where a
              * crossover should occur given the domain lengths, then the
              * junction must have one less stack.
              */
            bool check_junction_constraint(
                    Domain& cd_j2,
                    Domain& cd_j3,
                    Domain& cd_j4,
                    Domain& cd_k1,
                    Domain& cd_k2);
    };

    class ThreeQuarterTurnDomain: public Domain {
        using Domain::Domain;
        public:

            // Constraint checkers
            bool check_twist_constraint(VectorThree ndr, Domain& cd_j);
            bool check_kink_constraint(VectorThree ndr, Domain& cd_j);
            bool check_junction_constraint(
                    Domain& cd_j2,
                    Domain& cd_j3,
                    Domain& cd_j4,
                    Domain& cd_k1,
                    Domain& cd_k2);
    };
}

#endif // DOMAIN_H
