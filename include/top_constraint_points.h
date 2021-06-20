// cb_constraint_points.h

#ifndef TOP_CONSTRAINT_POINTS_H
#define TOP_CONSTRAINT_POINTS_H

#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "domain.h"
#include "ideal_random_walk.h"
#include "origami_system.h"
#include "parser.h"
#include "random_gens.h"
#include "utility.h"

namespace topConstraintPoints {

using std::pair;
using std::set;
using std::unordered_map;
using std::vector;

using domainContainer::Domain;
using idealRandomWalk::IdealRandomWalks;
using origami::OrigamiSystem;
using utility::VectorThree;

bool staple_excluded(vector<int> exclude_staples, int staple);
bool staple_excluded(set<int> exclude_staples, int staple);

/** Check if domain in given domains */
bool domain_included(vector<Domain*> domains, Domain* d);

/** Check if chain in given chain */
bool chain_included(vector<int> staples, int staple);
bool chain_included(set<int> staples, int staple);

/** Network of bound staples with potential growth- and endpoint info */
class StapleNetwork {
  public:
    StapleNetwork(OrigamiSystem& origami);

    void set_excluded_staples(vector<int> excluded_staples);

    /** Set internal scaffold domains */
    void set_scaffold_domains(vector<Domain*> scaffold_domains);

    /** Scan the staple network starting from the given staple domain
     *
     * This methods assumes that the given domain is bound to a
     * scaffold domain. The given scaffold segment is the regrowth
     * segment being considered, and is involved in determining
     * whether the network is considered external.
     */
    void scan_network(Domain* d);

    set<int> get_participating_chains();
    vector<pair<Domain*, Domain*>> get_potential_growthpoints();
    vector<pair<Domain*, Domain*>> get_potential_inactive_endpoints();
    vector<Domain*> get_potential_domain_stack();
    unordered_map<Domain*, int> get_staple_to_segs_map();
    unordered_map<pair<int, int>, int> get_dirs();
    bool externally_bound();

  private:
    void clear_network();
    void scan_staple_topology(Domain* domain);
    vector<Domain*> make_staple_stack(Domain* d, int ci);
    void add_potential_inactive_endpoint(Domain* d, Domain* bd);

    OrigamiSystem& m_origami;

    bool m_external {false};
    vector<Domain*> m_scaffold_ds; // Scaffold segment being regrown
    vector<int> m_ex_staples; // Excluded staples
    set<int> m_net_cs; // Chains involved in network
    vector<pair<Domain*, Domain*>> m_pot_gps {}; // Potential growthpoints
    vector<Domain*> m_pot_ds {}; // Potential domain stack
    vector<pair<Domain*, Domain*>> m_pot_iaes; // Potential inactive endpoints
    unordered_map<Domain*, int> m_segs {}; // Map from domain to segment

    // Map from domain and segment to regrowth direction
    unordered_map<pair<int, int>, int> m_domain_to_dir {};
};

/** Topology of an associated scaffold segment
 *
 * This class is intended to be reused for different scaffold segments. It
 * will store information on the binding states of the domains of given
 * scaffold segments and associated staples. It differentiates between
 * "externally bound" and "internally bound" staples. The former are those
 * that are also bound to a scaffold domain outside the set provided
 * (whether directly, or indirectly through other bound staples), while
 * the latter are the compliment. Growthpoints are associated with
 * internally bound staples, while externally bound staples are only
 * only associated with endpoints. A given domain of a chain can have
 * several endpoints associated with it. Active endpoints are those for
 * which one of the domains has a set configuration, while inactive are
 * those where neither domain has been set. Once a domain with an
 * inactive endpoint is set, it is converted to an active endpoint with
 * the current position (by calling update).
 */
class Constraintpoints {
  public:
    Constraintpoints(
            OrigamiSystem& origami_system,
            IdealRandomWalks& ideal_random_walks);

    /** Get direction of growth in current segment of current chain */
    int get_dir(Domain* d);

    /** Return endpoint positions erased from most recent update */
    vector<VectorThree> get_erased_endpoints();

    void reset_internal();

    /** Find all growthpoints and endpoints for given scaffold segment
     *
     * These include the active endpoints from the externally bound
     * staples and the inactive endpoints from the internally bound
     * staples. Note this does not include the endpoint of the
     * external scaffold domain that the scaffold domain is intended
     * to regrow to; the direction can be selected after this has been
     * called.
     */
    void calculate_constraintpoints(
            vector<Domain*> scaffold_domains,
            int dir,
            vector<int> excluded_staples);
    void calculate_constraintpoints(
            vector<vector<Domain*>> scaffold_segments,
            vector<int> dirs,
            vector<int> excluded_staples);

    /** Return the set of interally bound staples
     *
     * As they have no external anchor, these staples must be regrown
     * when the scaffold segment is regrown.
     */
    set<int> staples_to_be_regrown();

    /** Return domain stack to be regrown
     *
     * Includes first domain of given scaffold range even thought this
     * is generally not regrown by the movetypes as it is used as a
     * reference by the first domain being regrown.
     */
    vector<Domain*> domains_to_be_regrown();

    /** Test if given domain is used as a growthpoint */
    bool is_growthpoint(Domain* domain);

    /** Test if given domain grown out from a growthpoint */
    bool is_stemdomain(Domain* domain);

    /** Add active endpoint on given domain at given position */
    void add_active_endpoint(Domain* domain, VectorThree endpoint_pos);
    void add_active_endpoint(Domain* domain, VectorThree endpoint_pos, int seg);
    void add_inactive_endpoint(Domain* d_i, Domain* d_j);
    void add_growthpoint(Domain* growthpoint, Domain* stemd);
    void add_stem_seg_pair(Domain* stemd, vector<int> seg_pair);

    /** Reset the set of active endpoints to the initial set
     *
     * The initial set is taken to be the set after the last call
     * of calculate_constraintpoints
     */
    void reset_active_endpoints();

    /** Remove all active endpoints on the given domain */
    void remove_active_endpoint(Domain* domain);

    /** Remove all endpoints activated on the given domain */
    void remove_activated_endpoint(Domain* domain);

    /** Remove reached active endpoints and activate inactive ones */
    void update_endpoints(Domain* domain);

    /** Get domain that is to grow from given domain */
    Domain* get_domain_to_grow(Domain* domain);

    /** Get growthpoint of give stem domain */
    Domain* get_growthpoint(Domain* domain);

    /** Active endpoint reached */
    bool endpoint_reached(Domain* domain, VectorThree pos);

    /** Calculate the product of the number of IRW
     *
     * For the given domain, calculate the product of the number of
     * ideal random walks (IRW)s for each active endpoint on the
     * current chain. The number of remaining steps is taken from
     * the number of domains until the endpoint is reached, but an
     * offset can be given for this number.
     */
    long double calc_num_walks_prod(
            Domain* domain,
            VectorThree pos,
            int dir,
            int offset = 0);

    /** Return whether there any IRWs */
    bool walks_remain(Domain* domain, VectorThree pos, int offset = 0);

    // DEBUG
    vector<pair<int, VectorThree>> get_active_endpoints(int c_i, int seg);
    Domain* get_inactive_endpoints(Domain* domain);

  private:
    bool walks_remain(
            pair<int, int> key,
            Domain* domain,
            VectorThree pos,
            int dir,
            int offset = 0);
    void find_growthpoints_endpoints(
            vector<Domain*> scaffold_domains,
            vector<int> excluded_staples,
            int seg);
    bool bound_to_self(Domain* d);
    bool staple_already_checked(Domain* domain, set<int> checked_staples);
    void add_growthpoints(
            vector<pair<Domain*, Domain*>> potential_growthpoints);
    void add_inactive_endpoints(vector<pair<Domain*, Domain*>> pot_iaes);
    void add_regrowth_staples(
            set<int> participating_chains,
            vector<int> excluded_staples);
    void add_domains_to_stack(vector<Domain*> potential_d_stack);
    void add_active_endpoints_on_scaffold(
            vector<pair<Domain*, Domain*>> potential_growthpoints,
            vector<pair<Domain*, Domain*>> potential_inactive_endpoints,
            int seg);
    void add_staple_to_segs_maps(unordered_map<Domain*, int> s_seg_map);
    void add_dirs(unordered_map<pair<int, int>, int> domain_to_dirs);
    int calc_remaining_steps(
            int endpoint_d_i,
            Domain* domain,
            int dir,
            int step_offset);

    OrigamiSystem& m_origami_system;
    IdealRandomWalks& m_ideal_random_walks;
    StapleNetwork m_staple_network;

    // State variables
    // Scaffold domains to be regrown
    vector<Domain*> m_scaffold_domains {};

    // Staples to be regrown
    set<int> m_regrowth_staples;

    // Stack of domains in order for regrowth
    vector<Domain*> m_d_stack;

    // Chain ids of checked staples
    set<int> m_checked_staples;
    vector<VectorThree> m_erased_endpoints;

    // Map from domain to segment
    unordered_map<Domain*, int> m_segs {};

    // Map from domain and segment to regrowth direction
    unordered_map<pair<int, int>, int> m_domain_to_dir {};

    // Map from growthpoint domain to domain to grow
    unordered_map<Domain*, Domain*> m_growthpoints {};

    // Map from stem domain to growthpoint domain
    unordered_map<Domain*, Domain*> m_stemdomains {};

    // Map from chain to its active endpoints (absolute positions)
    unordered_map<pair<int, int>, vector<pair<int, VectorThree>>>
            m_active_endpoints {};

    // Map from chain to its inactive endpoints
    // (domain that hasn't been set)
    unordered_map<pair<int, int>, vector<pair<int, VectorThree>>>
            m_initial_active_endpoints {};

    // Map from chain to endpoints it will impose once grown
    unordered_map<Domain*, Domain*> m_inactive_endpoints {};

    // Stem domain to segments to grow out
    unordered_map<Domain*, vector<int>> m_stemd_to_segs {};
};
} // namespace topConstraintPoints

#endif // TOP_CONSTRAINT_POINTS
