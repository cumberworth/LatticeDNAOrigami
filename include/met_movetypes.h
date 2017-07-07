// met_movetypes.h

#ifndef MET_STAPLE_MOVETYPES_H
#define MET_STAPLE_MOVETYPES_H

#include "movetypes.h"


namespace MetMovetypes {

    using namespace Movetypes;

    class MetMCMovetype: public RegrowthMCMovetype {
        public:
            using RegrowthMCMovetype::RegrowthMCMovetype;

            void reset_internal();
            string m_label() {return "MetMCMovetype";};
        protected:
            double m_delta_e {0};

            void grow_chain(vector<Domain*> domains);
            void unassign_domains(vector<Domain*> domains);

            void update_bias(int sign);
    };

    class MetStapleExchangeMCMovetype: public MetMCMovetype {
        public:
            using MetMCMovetype::MetMCMovetype;
            bool attempt_move();

            string m_label() {return "MetStapleExchangeMCMovetype";}
            void reset_internal();
        protected:

            // These can be overidden for a derived class the excludes misbinding
            int preconstrained_df {1};

            // I select by domains not by sites, so this is consistent
            int m_insertion_sites {m_origami_system.num_domains()};

        private:
            bool staple_insertion_accepted(int c_i_ident);
            bool staple_deletion_accepted(int c_i_ident);
            bool insert_staple();
            bool delete_staple();
    };

    class MetStapleRegrowthMCMovetype: public MetMCMovetype {
        public:
            using MetMCMovetype::MetMCMovetype;
            bool attempt_move();

            string m_label() {return "MetStapleRegrowthMCMovetype";}
    };
}

#endif // MET_MOVETYPES_H
