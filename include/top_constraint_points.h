// cb_constraint_points.h

#ifndef TOP_CONSTRAINT_POINTS_H
#define TOP_CONSTRAINT_POINTS_H

#include <vector>
#include <set>
#include <unordered_map>
#include <utility>

#include "ideal_random_walk.h"
#include "domain.h"
#include "origami_system.h"
#include "parser.h"
#include "random_gens.h"
#include "utility.h"

namespace TopConstraintPoints {

    using std::pair;
    using std::set;
    using std::unordered_map;
    using std::vector;

    using DomainContainer::Domain;
    using IdealRandomWalk::IdealRandomWalks;
    using Origami::OrigamiSystem;
    using Parser::InputParameters;
    using RandomGen::RandomGens;
    using Utility::VectorThree;

    class Constraintpoints {
        public:
            Constraintpoints(
                    OrigamiSystem& origami_system,
                    IdealRandomWalks& ideal_random_walks);

            void reset_internal();
            void calculate_constraintpoints(vector<Domain*> scaffold_domains);
            inline set<int> staples_to_be_regrown() {return m_regrowth_staples;}
            bool is_growthpoint(Domain* domain);
            void add_active_endpoint(Domain* domain, VectorThree endpoint_pos);
            void reset_active_endpoints();
            void remove_active_endpoint(Domain* domain);
            void update_endpoints(Domain* domain);
            Domain* get_domain_to_grow(Domain* domain);
            bool endpoint_reached(Domain* domain, VectorThree pos);
            long double calc_num_walks_prod(Domain* domain, VectorThree pos,
                    vector<Domain*> domains, int dir, int step_offset=0);

            // For debugging
            inline vector<pair<int, VectorThree>> get_active_endpoints(int c_i) {
                    return m_active_endpoints[c_i];}
            inline Domain* get_inactive_endpoints(Domain* domain) {
                return m_inactive_endpoints[domain];}

        private:
            void find_staples_growthpoints_endpoints(vector<Domain*> scaffold_domains);
            bool staple_already_checked(Domain* domain, set<int> checked_staples);
            void add_growthpoints(vector<pair<Domain*, Domain*>> potential_growthpoints);
            void add_regrowth_staples(set<int> participating_chains);
            void add_active_endpoints_on_scaffold(
                    set<int> participating_chains,
                    vector<pair<Domain*, Domain*>> potential_growthpoints);
            void scan_staple_topology(
                    Domain* domain,
                    set<int>& potential_growthpoints,
                    vector<pair<Domain*, Domain*>>& participating_chains,
                    vector<Domain*>& scaffold_domains,
                    bool& externally_bound);
            int calc_remaining_steps(int endpoint_d_i, Domain* domain, int dir,
                    int step_offset);

            OrigamiSystem& m_origami_system;
            IdealRandomWalks& m_ideal_random_walks;

            // Scaffold domains to be regrown
            vector<Domain*> m_scaffold_domains {};

            // Staples to be regrown
            set<int> m_regrowth_staples;

            // Map from growthpoint domain to domain to grow
            unordered_map<Domain*, Domain*> m_growthpoints {};

            // Map from chain to it's active endpoints
            unordered_map<int, vector<pair<int, VectorThree>>> m_active_endpoints {};
            unordered_map<int, vector<pair<int, VectorThree>>> m_initial_active_endpoints {};

            // Map from chain to endpoints it will impose once grown
            unordered_map<Domain*, Domain*> m_inactive_endpoints {};

    };
}

#endif // TOP_CONSTRAINT_POINTS
