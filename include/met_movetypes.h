// met_movetypes.h

#ifndef MET_STAPLE_MOVETYPES_H
#define MET_STAPLE_MOVETYPES_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "domain.h"
#include "hash.h"
#include "movetypes.h"
#include "utility.h"

namespace movetypes {

    using std::ostream;
    using std::string;
    using std::unordered_map;
    using std::vector;

    using domainContainer::Domain;
    using movetypes::MovetypeTracking;
    using movetypes::RegrowthMCMovetype;
    using utility::StapleExchangeTracking;
    using utility::StapleRegrowthTracking;

    class MetMCMovetype:
        virtual public RegrowthMCMovetype {

        public:
            MetMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params);
            MetMCMovetype(const MetMCMovetype&) = delete;
            MetMCMovetype& operator=(const MetMCMovetype&) = delete;
            void reset_internal() override;

        protected:
            void grow_chain(vector<Domain*> domains) override;
            void add_external_bias() override;

            void unassign_domains(vector<Domain*> domains);

            double m_delta_e {0};
    };

    class MetStapleExchangeMCMovetype:
        public MetMCMovetype {

        public:
            MetStapleExchangeMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params,
                    vector<double> exchange_mults,
                    bool adaptive_exchange);
            MetStapleExchangeMCMovetype(const
                    MetStapleExchangeMCMovetype&) = delete;
            MetStapleExchangeMCMovetype& operator=(const
                    MetStapleExchangeMCMovetype&) = delete;

            void reset_internal() override;
            void write_log_summary(ostream* log_stream) override;

        private:
            bool internal_attempt_move() override;
            void add_tracker(bool accepted) override;

            bool staple_insertion_accepted(int c_i_ident);
            bool staple_deletion_accepted(int c_i_ident);
            bool insert_staple();
            bool delete_staple();

            bool m_adaptive_exchange;
            bool m_staple_bound {false};
            vector<double> m_exchange_mults;
            StapleExchangeTracking m_tracker {};
            unordered_map<StapleExchangeTracking, MovetypeTracking> m_tracking {};

            // These need to be changed if misbinding is excluded
            int preconstrained_df {1};

            // I select by domains not by sites, so this is consistent
            int m_insertion_sites {m_origami_system.num_domains()};

    };

    class MetStapleRegrowthMCMovetype:
        public MetMCMovetype {

        public:
            MetStapleRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params);
            MetStapleRegrowthMCMovetype(const
                    MetStapleRegrowthMCMovetype&) = delete;
            MetStapleRegrowthMCMovetype& operator=(const
                    MetStapleRegrowthMCMovetype&) = delete;

            void write_log_summary(ostream* log_stream) override;

        private:
            bool internal_attempt_move() override;
            void add_tracker(bool accepted) override;

            StapleRegrowthTracking m_tracker {};
            unordered_map<StapleRegrowthTracking, MovetypeTracking> m_tracking {};
    };
}

#endif // MET_MOVETYPES_H
