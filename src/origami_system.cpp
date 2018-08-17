// origami_system.cpp

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>

#include "bias_functions.h"
#include "order_params.h"
#include "origami_system.h"
#include "files.h"

namespace origami {

    using std::abs;
    using std::cout;
    using std::max_element;

    using biasFunctions::SystemBiases;
    using files::OrigamiInputFile;
    using files::OrigamiTrajInputFile;
    using domainContainer::SixteenDomain;
    using orderParams::SystemOrderParams;
    using utility::NotImplemented;
    using utility::OrigamiMisuse;

    bool Chain::operator==(Chain chain_2) {
        bool index_match {this->index == chain_2.index};
        bool identity_match {this->identity == chain_2.index};
        bool positions_match {this->positions == chain_2.positions};
        bool orientations_match {this->orientations == chain_2.orientations};
        if (index_match and identity_match and positions_match and orientations_match) {
            return true;
        }
        else {
            return false;
        }
    }

    OrigamiSystem::OrigamiSystem(
            const vector<vector<int>>& identities,
            const vector<vector<string>>& sequences,
            const Chains& chains,
            bool cyclic,
            double volume,
            double staple_u,
            InputParameters& params) :

            m_identities {identities},
            m_sequences {sequences},
            m_temp {params.m_temp},
            m_volume {volume},
            m_cation_M {params.m_cation_M},
            m_staple_u {staple_u},
            m_cyclic {cyclic},
            m_pot {identities, sequences, params} {

        initialize_complementary_associations();
        initialize_scaffold(chains[0]);
        initialize_staples(chains);
        set_all_domains(chains);

        m_ops = std::make_unique<SystemOrderParams>(params, *this);
        m_biases = std::make_unique<SystemBiases>(*this, *m_ops, params);
    }

    OrigamiSystem::~OrigamiSystem() {
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                delete domain;
            }
        }
    }

    vector<vector<Domain*>> OrigamiSystem::get_chains() {
        return m_domains;
    }

    vector<Domain*> OrigamiSystem::get_last_chain() {
        return m_domains.back();
    }

    int OrigamiSystem::num_staples() const {
        return m_num_staples;
    }

    int OrigamiSystem::num_domains() {
        return m_num_domains;
    }

    int OrigamiSystem::num_bound_domain_pairs() const {
        return m_num_bound_domain_pairs;
    }

    int OrigamiSystem::num_fully_bound_domain_pairs() const {
        return m_num_fully_bound_domain_pairs;
    }

    int OrigamiSystem::num_self_bound_domain_pairs() const {
        return m_num_self_bound_domain_pairs;
    }

    int OrigamiSystem::num_misbound_domain_pairs() const {
        return num_bound_domain_pairs() - num_fully_bound_domain_pairs();
    }

    int OrigamiSystem::num_stacked_domain_pairs() const {
        return m_num_stacked_domain_pairs;
    }

    int OrigamiSystem::num_linear_helix_trips() const {
        return m_num_linear_helix_trips;
    }

    int OrigamiSystem::num_stacked_junct_quads() const {
        return m_num_stacked_junct_quads;
    }

    int OrigamiSystem::num_staples_of_ident(int staple_ident) const {
        return m_identity_to_index[staple_ident].size();
    }

    vector<int> OrigamiSystem::staples_of_ident(int c_ident) {
        return m_identity_to_index[c_ident];
    }

    vector<int> OrigamiSystem::complementary_scaffold_domains(int staple_ident)
            const {
        return m_staple_ident_to_scaffold_ds[staple_ident];
    }

    Domain* OrigamiSystem::unbound_domain_at(VectorThree pos) const {
        return m_pos_to_unbound_d.at(pos);
    }

    bool OrigamiSystem::check_domains_complementary(Domain& cd_i, Domain& cd_j) {
        return m_pot.check_domains_complementary(cd_i, cd_j);
    }

    double OrigamiSystem::energy() const {
        return m_energy;
    }

    SystemOrderParams& OrigamiSystem::get_system_order_params() {
        return *m_ops;
    }

    SystemBiases& OrigamiSystem::get_system_biases() {
        return *m_biases;
    }

    vector<Domain*> OrigamiSystem::get_chain(int c_i) {
        int c_i_index {utility::index(m_chain_indices, c_i)};
        return m_domains[c_i_index];
    }

    Domain* OrigamiSystem::get_domain(int c_i, int d_i) {
        return get_chain(c_i)[d_i];
    }

    vector<int> OrigamiSystem::get_staple_counts() {
        vector<int> staple_counts(m_identities.size() - 1, 0);
        for (size_t c_i {1}; c_i != m_domains.size(); c_i++) {
            int c_ident {m_domains[c_i][0]->m_c_ident};
            staple_counts[c_ident - 1]++;
        }

        return staple_counts;
    }

    int OrigamiSystem::num_unique_staples() const {
        int unique_staple_count {0};
        for (auto indices: m_identity_to_index) {
            if (indices.size() > 0) {
                unique_staple_count++;
            }
        }

        // Minus one to remove scaffold
        return unique_staple_count - 1;
    }

    Chains OrigamiSystem::chains() const {
        // Return chains data structure for current config
        Chains chains;
        for (size_t i {0}; i != m_domains.size(); i++) {
            int c_i_ident {m_chain_identities.at(i)};
            int c_i {m_chain_indices.at(i)};
            vector<VectorThree> positions;
            vector<VectorThree> orientations;
            for (auto domain: m_domains[i]) {
                positions.push_back(domain->m_pos);
                orientations.push_back(domain->m_ore);
            }
            Chain chain {c_i, c_i_ident, positions, orientations};
            chains.push_back(chain);
        }
        return chains;
    }

    Occupancy OrigamiSystem::position_occupancy(VectorThree pos) const {
        Occupancy occ;
        if (m_position_occupancies.count(pos) == 0) {
            occ = Occupancy::unassigned;
        }
        else {
            occ = m_position_occupancies.at(pos);
        }
        return occ;
    }

    void OrigamiSystem::update_enthalpy_and_entropy() {
        m_hyb_enthalpy = 0;
        m_hyb_entropy = 0;
        set<Domain*> accounted_domains; // Domains whose energy is already counted
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                Occupancy state {domain->m_state};
                if (state == Occupancy::bound or state == Occupancy::misbound) {
                    if (find(accounted_domains.begin(), accounted_domains.end(),
                                domain) == accounted_domains.end()) {
                        Domain& bound_domain {*domain->m_bound_domain};

                        m_hyb_enthalpy += m_pot.hybridization_enthalpy(
                                *domain, bound_domain);
                        m_hyb_entropy += m_pot.hybridization_entropy(
                                *domain, bound_domain);
                        accounted_domains.insert(domain->m_bound_domain);
                    }
                }
            }
        }
        m_stacking_energy = m_energy - (m_hyb_enthalpy - m_hyb_entropy);
    }

    double OrigamiSystem::hybridization_enthalpy() {
        return m_hyb_enthalpy;
    }

    double OrigamiSystem::hybridization_entropy() {
        return m_hyb_entropy;
    }

    double OrigamiSystem::stacking_energy() {
        return m_stacking_energy;
    }

    bool OrigamiSystem::configuration_fully_set() {
        if (m_num_unassigned_domains == 0) {
            return true;
        }
        else {
            return false;
        }
    }

    int OrigamiSystem::num_unassigned_domains() {
        return m_num_unassigned_domains;
    }

    void OrigamiSystem::check_all_constraints() {

        // Unassign everything (and check nothing was already unassigned)
        if (m_num_unassigned_domains != 0) {
            throw OrigamiMisuse {};
        }
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                if (domain->m_state == Occupancy::unassigned) {
                    cout << "Domain unassigned after move complete\n";
                    throw OrigamiMisuse {};
                }
                else {
                    unassign_domain(*domain);
                }
            }
        }

        // Check that number of stacked pairs is consistent
        if (m_num_stacked_domain_pairs != 0) {
            cout << "Stack count inconsistency\n";
            cout << m_num_stacked_domain_pairs << "\n";
            set_all_domains();
            throw OrigamiMisuse {};
        }

        // Check that energy has returned to 0 (eps is totally arbitrary)
        double eps {0.000001};
        if (m_energy < -eps or m_energy > eps) {
            cout << "Inconsistency in system energy\n";
            cout << m_energy << "\n";
            set_all_domains();
            throw OrigamiMisuse {};
        }
        else {

            // Prevent rounding errors from building up
            m_energy = 0;
        }

        // Reset configuration
        set_all_domains();
    }

    void OrigamiSystem::set_config(Chains chains) {

        // Unassign and delete everything but scaffold
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                unassign_domain(*domain);
            }
        }
        for (int i {static_cast<int>(m_domains.size() - 1)}; i != 0; i--) {
            delete_chain(m_domains[i][0]->m_c);
        }

        initialize_staples(chains);
        set_all_domains(chains);
    }

    double OrigamiSystem::check_domain_constraints(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore) {

        DeltaConfig delta_config {internal_check_domain_constraints(
                cd_i, pos, ore)};
        internal_unassign_domain(cd_i);

        return delta_config.e;
    }

    void OrigamiSystem::check_distance_constraints() {
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                Domain* next_domain {*domain + 1};
                if (next_domain == nullptr) {
                    continue;
                }
                else {
                    VectorThree dist {next_domain->m_pos - domain->m_pos};
                    if (dist.abssum() != 1) {
                        cout << "Contiguous domains not on adjacent sites\n";
                        throw OrigamiMisuse {};
                    }
                }
            }
        }
    }

    double OrigamiSystem::unassign_domain(Domain& cd_i) {
        DeltaConfig delta_config {internal_unassign_domain(cd_i)};
        m_energy += delta_config.e;
        m_num_stacked_domain_pairs += (delta_config.stacked_pairs/2);
        m_num_linear_helix_trips += delta_config.linear_helices;
        m_num_stacked_junct_quads += delta_config.stacked_juncts;
        m_num_unassigned_domains++;

        return delta_config.e;
    }

    int OrigamiSystem::add_chain(int c_i_ident) {
        // Add chain with domains in unassigned state and return chain index.
        m_current_c_i += 1;
        return add_chain(c_i_ident, m_current_c_i);
    }

    int OrigamiSystem::add_chain(int c_i_ident, int c_i) {
        // Add chain with given index
        m_identity_to_index[c_i_ident].push_back(c_i);
        m_chain_indices.push_back(c_i);
        m_chain_identities.push_back(c_i_ident);

        int chain_length {static_cast<int>(m_identities[c_i_ident].size())};
        m_domains.push_back({});
        m_num_staples++;
        Domain* prev_domain {nullptr};
        for (int d_i {0}; d_i != chain_length; d_i++) {
            int d_i_ident {m_identities[c_i_ident][d_i]};
            Domain* domain;
            // Assume that all domains are 16
//            size_t domain_size {m_sequences[c_i_ident][d_i].size()};
//            if (domain_size == 16 or domain_size == 15) {
            domain = new SixteenDomain {c_i, c_i_ident, d_i, d_i_ident, chain_length};
//            }
 //           else {
  //              throw NotImplemented {};
   //         }

            // Set forward and backwards domains
            if (prev_domain != nullptr) {
                prev_domain->m_forward_domain = domain;
            }

            domain->m_backward_domain = prev_domain;
            domain->m_forward_domain = nullptr;

            m_domains.back().push_back(domain);
            prev_domain = domain;

            m_num_domains++;
            m_num_unassigned_domains++;
        }

        return c_i;
    }

    void OrigamiSystem::delete_chain(int c_i) {
        // Delete chain c_i_
        
        int c_i_index {utility::index(m_chain_indices, c_i)};
        int c_i_ident {m_chain_identities[c_i_index]};

        // Index in m_identity_to_index of given index and type
        int j {utility::index(m_identity_to_index[c_i_ident], c_i)};
        m_identity_to_index[c_i_ident].erase(m_identity_to_index[c_i_ident].begin()
                + j);
        m_chain_indices.erase(m_chain_indices.begin() + c_i_index);
        m_chain_identities.erase(m_chain_identities.begin() + c_i_index);
        m_num_domains -= m_domains[c_i_index].size();
        m_num_staples--;
        for (auto domain: m_domains[c_i_index]) {
            delete domain;
            m_num_unassigned_domains--;
        }
        m_domains.erase(m_domains.begin() + c_i_index);
    }

    void OrigamiSystem::temp_reduce_staples_by_one() {
        m_num_staples--;
    }

    void OrigamiSystem::undo_reduce_staples_by_one() {
        m_num_staples++;
    }

    double OrigamiSystem::set_checked_domain_config(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore) {
        update_domain(cd_i, pos, ore);
        update_occupancies(cd_i, pos);

        // Update internal energy
        double delta_e {0};
        if (cd_i.m_state == Occupancy::misbound) {
            delta_e += m_pot.hybridization_energy(cd_i, *cd_i.m_bound_domain);
        }
        else if (cd_i.m_state == Occupancy::bound) {
            delta_e += m_pot.hybridization_energy(cd_i, *cd_i.m_bound_domain);
            DeltaConfig delta_config;
            delta_config = m_pot.check_stacking(cd_i, *cd_i.m_bound_domain);
            delta_e += delta_config.e;
            m_num_stacked_domain_pairs += (delta_config.stacked_pairs/2);
            m_num_linear_helix_trips += delta_config.linear_helices;
            m_num_stacked_junct_quads += delta_config.stacked_juncts;
        }
        m_energy += delta_e;
        m_num_unassigned_domains--;
        return delta_e;
    }

    double OrigamiSystem::set_domain_config(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore) {
        if (cd_i.m_state != Occupancy::unassigned) {
            cout << "Trying to set an already assigned domain\n";
            throw OrigamiMisuse {};
        }

        // Check constraints and update if obeyed, otherwise throw
        DeltaConfig delta_config {internal_check_domain_constraints(
                cd_i, pos, ore)};
        if (m_constraints_violated) {
            internal_unassign_domain(cd_i);
        }
        else{
            m_energy += delta_config.e;
            m_num_stacked_domain_pairs += (delta_config.stacked_pairs/2);
            m_num_linear_helix_trips += delta_config.linear_helices;
            m_num_stacked_junct_quads += delta_config.stacked_juncts;
            m_num_unassigned_domains--;
        }
        return delta_config.e;
    }

    void OrigamiSystem::set_domain_orientation(Domain& cd_i, VectorThree ore) {
        Occupancy occupancy {cd_i.m_state};
        if (occupancy == Occupancy::bound) {
            m_constraints_violated = true;
        }
        else {
            cd_i.m_ore = ore;
        }
    }

    void OrigamiSystem::center(int centering_domain) {
        // Translate the system such that the first scaffold domain is on the origin

        // The maps that go from a position to something need to be reset
        unordered_map<VectorThree, Domain*> pos_to_unbound_d {};
        unordered_map<VectorThree, Occupancy> position_occupancies {};
        VectorThree refpos {m_domains[0][centering_domain]->m_pos};
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                VectorThree new_pos {domain->m_pos - refpos};
                pos_to_unbound_d[new_pos] = m_pos_to_unbound_d[domain->m_pos];
                position_occupancies[new_pos] = m_position_occupancies[domain->m_pos];
                domain->m_pos = new_pos;
            }
        }
        m_pos_to_unbound_d = pos_to_unbound_d;
        m_position_occupancies = position_occupancies;
    }

    void OrigamiSystem::set_all_domains() {
        // Use current positions and orientations
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                set_domain_config(*domain, domain->m_pos, domain->m_ore);
                if (m_constraints_violated) {
                    cout << "Constaints in violation after move complete\n";
                    set_domain_config(*domain, domain->m_pos, domain->m_ore);
                    throw OrigamiMisuse {};
                }
            }
        }
        check_distance_constraints();
    }

    void OrigamiSystem::set_all_domains(Chains config) {
        // Use given positions and orientations

        // Set position and orientation of domains
        for (size_t i {0}; i != config.size(); i++) {
            Chain chain {config[i]};
            int num_domains {static_cast<int>(m_domains[i].size())};
            for (int d_i {0}; d_i != num_domains; d_i++) {
                Domain* domain {m_domains[i][d_i]};
                VectorThree pos = chain.positions[d_i];
                VectorThree ore = chain.orientations[d_i];
                domain->m_pos = pos;
                domain->m_ore = ore;
            }
        }

        // Set all domains and check all constraints
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                set_domain_config(*domain, domain->m_pos, domain->m_ore);
                if (m_constraints_violated) {
                    cout << "Constaints in violation after move complete\n";
                    set_domain_config(*domain, domain->m_pos, domain->m_ore);
                    throw OrigamiMisuse {};
                }
            }
        }
        check_distance_constraints();
    }

    void OrigamiSystem::update_temp(double temp) {
        m_temp = temp;
        
        m_pot.update_temp(temp);

        // Recalculate system energy
        update_energy();

        // Update volume such that it cancels fugacity
        m_volume = exp(-m_staple_u / temp);
    }

    void OrigamiSystem::update_staple_u(double u) {
        m_staple_u = u;

        // Update volume such that it cancels fugacity
        m_volume = exp(-u / m_temp);
    }

    void OrigamiSystem::initialize_complementary_associations() {
        // Intialize staple identity to complementary scaffold domains container
        // Staple identities are 1 indexed (scaffold is 0)
        m_staple_ident_to_scaffold_ds.push_back({});
        m_identity_to_index.push_back({});
        for (unsigned int i {1}; i != m_identities.size(); ++i) {
            m_identity_to_index.push_back({});
            vector<int> staple {m_identities[i]};
            vector<int> scaffold_d_is {};
            for (auto staple_d_i: staple) {

                // staple_d_i's are negatives of scaffold_d_i's
                if (std::find(m_identities[c_scaffold].begin(),
                        m_identities[c_scaffold].end(), -staple_d_i) ==
                        m_identities[c_scaffold].end()) {

                    scaffold_d_is.push_back(0);
                }
                else {
                    int scaffold_d_i {utility::index(m_identities[c_scaffold],
                            -staple_d_i)};
                    scaffold_d_is.push_back(scaffold_d_i);
                }
            }
            m_staple_ident_to_scaffold_ds.push_back(scaffold_d_is);
        }
    }

    void OrigamiSystem::initialize_staples(Chains chains) {

        // Create domain objects
        for (size_t i {1}; i != chains.size(); i++) {
            Chain chain {chains[i]};
            int c_i {chain.index};
            int c_i_ident {chain.identity};
            add_chain(c_i_ident, c_i);
        }

        // Current unique chain index
        m_current_c_i = *max_element(m_chain_indices.begin(),
                m_chain_indices.end());
    }

    void OrigamiSystem::initialize_scaffold(Chain scaffold_chain) {
        int c_i {scaffold_chain.index};
        int c_i_ident {scaffold_chain.identity};
        add_chain(c_i_ident, c_i);
        m_num_staples--;

        // Make scaffold chain domains modular if cyclic
        if (m_cyclic) {
            Domain* first_domain {m_domains[c_i][0]};
            Domain* last_domain {m_domains[c_i].back()};
            last_domain->m_forward_domain = first_domain;
            first_domain->m_backward_domain = last_domain;
        }
    }

    DeltaConfig OrigamiSystem::internal_unassign_domain(Domain& cd_i) {
        // Deletes positions, orientations, and removes/unassigns occupancies.
        Occupancy occupancy {cd_i.m_state};
        DeltaConfig delta_config {};
        switch (occupancy) {
            case Occupancy::bound:
                m_num_fully_bound_domain_pairs -= 1;
                m_num_bound_domain_pairs -= 1;
                delta_config = m_pot.check_stacking(cd_i, *cd_i.m_bound_domain);
                delta_config.e *= -1;
                delta_config.stacked_pairs *= -1;
                delta_config.linear_helices *= -1;
                delta_config.stacked_juncts *= -1;
                delta_config.e += unassign_bound_domain(cd_i);
                break;
            case Occupancy::misbound:
                m_num_bound_domain_pairs -= 1;
                if (cd_i.m_bound_domain->m_c == cd_i.m_c) {
                    m_num_self_bound_domain_pairs -= 1;
                }
                delta_config.e += unassign_bound_domain(cd_i);
                break;
            case Occupancy::unbound:
                unassign_unbound_domain(cd_i);
                break;
            case Occupancy::unassigned:

                // Allows for double unassignment
                m_num_unassigned_domains--;
                break;
        }
        return delta_config;
    }


    double OrigamiSystem::unassign_bound_domain(Domain& cd_i) {
        Domain& cd_j {*cd_i.m_bound_domain};
        double delta_e {-m_pot.hybridization_energy(cd_i, cd_j)};

        cd_i.m_bound_domain = nullptr;
        cd_j.m_bound_domain = nullptr;
        cd_i.m_state = Occupancy::unassigned;

        m_pos_to_unbound_d[cd_i.m_pos] = &cd_j;
        m_position_occupancies[cd_i.m_pos] = Occupancy::unbound;
        cd_j.m_state = Occupancy::unbound;
        return delta_e;
    }

    void OrigamiSystem::unassign_unbound_domain(Domain& cd_i) {
        VectorThree pos {cd_i.m_pos};
        m_pos_to_unbound_d.erase(pos);
        m_position_occupancies.erase(pos);
        cd_i.m_state = Occupancy::unassigned;
    }

    void OrigamiSystem::update_domain(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore) {
        cd_i.m_pos = pos;
        cd_i.m_ore = ore;
    }    

    void OrigamiSystem::update_occupancies(Domain& cd_i, VectorThree pos) {
        Occupancy occupancy {position_occupancy(pos)};
        Occupancy new_state;

        switch (occupancy) {
            case Occupancy::unbound: {
                Domain* cd_j {unbound_domain_at(pos)};
                m_num_bound_domain_pairs += 1;
                if (cd_i.m_d_ident == -cd_j->m_d_ident) {
                    new_state = Occupancy::bound;
                    m_num_fully_bound_domain_pairs += 1;
                }
                else {
                    if (cd_i.m_c == cd_j->m_c) {
                        m_num_self_bound_domain_pairs += 1;
                    }
                    new_state = Occupancy::misbound;
                }

                m_pos_to_unbound_d.erase(pos);
                cd_i.m_state = new_state;
                cd_j->m_state = new_state;
                m_position_occupancies[pos] = new_state;
                cd_j->m_bound_domain = &cd_i;
                cd_i.m_bound_domain = cd_j;
                break;
            }
            case Occupancy::unassigned:
                new_state = Occupancy::unbound;
                cd_i.m_state = new_state;
                m_position_occupancies[pos] = new_state;
                m_pos_to_unbound_d[pos] = &cd_i;
                break;
            default:
                cout << "Trying to bind to an already bound domain\n";
                throw OrigamiMisuse {};
        }
    }

    void OrigamiSystem::update_energy() {
        // Recalculate total system energy with new energy tables
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                unassign_domain(*domain);
            }
        }
        m_energy = 0;
        set_all_domains();
    }

    DeltaConfig OrigamiSystem::internal_check_domain_constraints(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore) {
        // Updates positions and orientations and returns without reverting if no 
        // constraint violation. But states left unassigned.
        Occupancy occupancy {position_occupancy(pos)};
        DeltaConfig delta_config {};

        switch (occupancy) {
            case Occupancy::bound:
                m_constraints_violated = true;
                m_num_unassigned_domains++;
                break;
            case Occupancy::misbound:
                m_constraints_violated = true;
                m_num_unassigned_domains++;
                break;
            case Occupancy::unbound:
                update_domain(cd_i, pos, ore);
                update_occupancies(cd_i, pos);
                delta_config = m_pot.bind_domain(cd_i);
                m_constraints_violated = m_pot.m_constraints_violated;
                break;
            case Occupancy::unassigned:
                update_domain(cd_i, pos, ore);
                update_occupancies(cd_i, pos);
        }
        return delta_config;
    }

    OrigamiSystemWithBias::OrigamiSystemWithBias(
            const vector<vector<int>>& identities,
            const vector<vector<string>>& sequences,
            const Chains& chains,
            bool cyclic,
            double volume,
            double staple_u,
            InputParameters& params) :
            OrigamiSystem(
                    identities,
                    sequences,
                    chains,
                    cyclic,
                    volume,
                    staple_u,
                    params) {
    }

    double OrigamiSystemWithBias::check_domain_constraints(Domain& cd_i,
            VectorThree pos, VectorThree ore) {

        double delta_e {OrigamiSystem::check_domain_constraints(cd_i, pos, ore)};

        // Determine what state would be if bound
        Occupancy state;
        Occupancy pos_state {position_occupancy(pos)};
        if (pos_state == Occupancy::unassigned) {
            state = Occupancy::unbound;
        }
        else if (pos_state == Occupancy::unbound) {
            if (m_pot.check_domains_complementary(*unbound_domain_at(pos),
                    cd_i)) {
                state = Occupancy::bound;
            }
            else {
                state = Occupancy::misbound;
            }
        }
        else {
            state = Occupancy::unassigned;
        }

        m_ops->check_one_domain(cd_i, pos, ore, state);
        delta_e += m_biases->check_one_domain(cd_i);

        return delta_e;
    }

    void OrigamiSystemWithBias::delete_chain(int c_i) {
        OrigamiSystem::delete_chain(c_i);
        // Potentially update ops here
    }

    double OrigamiSystemWithBias::unassign_domain(Domain& cd_i) {
        double delta_e {OrigamiSystem::unassign_domain(cd_i)};
        m_ops->update_one_domain(cd_i);
        delta_e += m_biases->calc_one_domain(cd_i);

        return delta_e;
    }

    double OrigamiSystemWithBias::set_checked_domain_config(Domain& cd_i,
            VectorThree pos, VectorThree ore) {
        double delta_e {OrigamiSystem::set_checked_domain_config(cd_i, pos, ore)};

        // If you change this to somehow use the checked value, lot's to check
        m_ops->update_one_domain(cd_i);
        delta_e += m_biases->calc_one_domain(cd_i);

        return delta_e;
    }

    double OrigamiSystemWithBias::set_domain_config(Domain& cd_i, VectorThree pos,
            VectorThree ore) {
        double delta_e {OrigamiSystem::set_domain_config(cd_i, pos, ore)};
        delta_e += m_biases->calc_one_domain(cd_i);

        return delta_e;
    }

    //void set_domain_orientation(Domain& cd_i, VectorThree ore) {
    //    delta_e {OrigamiSystem::set_domain_config(cd_i, pos, ore)};
    //    m_ops->update_one_domain(cd_i);
    //    delta_e += m_biases->calc_one_domain(cd_i);
    //}

    double molarity_to_lattice_volume(double molarity) {
        // Volume is in units of number of lattice sites.
        // u = KB*T*ln(p), where p is the number of particles per lattice site
        // z = exp(1/(KB*T)*u) = exp(ln(p)) = p
        // V * p = 1, V = 1 / p
        // So just convert molarity to number of particles per lattice site
        // The lattice site volume removes the units
        double V {1 / molarity};
        return V;
    }

    double molarity_to_chempot(double molarity, double temp) {
        double chempot {temp * log(molarity)};
        return chempot;
    }

    double chempot_to_volume(double chempot, double temp) {
        return exp(-chempot / temp);
    }

    OrigamiSystem* setup_origami(InputParameters& params) {

        OrigamiInputFile origami_input {params.m_origami_input_filename};
        vector<vector<int>> identities {origami_input.get_identities()};
        vector<vector<string>> sequences {origami_input.get_sequences()};
        vector<Chain> configs = origami_input.get_config();
        bool cyclic {origami_input.is_cyclic()};
        if (params.m_restart_traj_file != "") {
            OrigamiTrajInputFile traj_file {params.m_restart_traj_file};
            configs = traj_file.read_config(params.m_restart_step);
        }

        // Calculate chemical potential from specified staple concentration
        double staple_u {molarity_to_chempot(params.m_staple_M,
                params.m_temp_for_staple_u)};
        staple_u *= params.m_staple_u_mult;
        double volume {chempot_to_volume(staple_u, params.m_temp)};

        OrigamiSystem* origami;
        if (params.m_domain_update_biases_present) {
            origami = new OrigamiSystemWithBias {identities, sequences, configs,
                    cyclic, volume, staple_u, params};
        }
        else {
            origami = new OrigamiSystem {identities, sequences, configs, cyclic,
                volume, staple_u, params};
        }

        return origami;
    }
}
