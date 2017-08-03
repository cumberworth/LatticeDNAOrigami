// origami_system.cpp

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/utility.hpp>

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
            m_energy_filebase {params.m_energy_filebase} {

        initialize_complementary_associations();
        initialize_scaffold(chains[0]);
        initialize_staples(chains);
        get_energies();
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

    bool OrigamiSystem::check_domains_complementary(Domain& cd_i, Domain& cd_j) {
        bool comp;
        if (cd_i.m_d_ident == -cd_j.m_d_ident) {
            comp = true;
        }
        else {
            comp = false;
        }
        return comp;
    }

    ThermoOfHybrid OrigamiSystem::enthalpy_and_entropy() {
        ThermoOfHybrid DH_DS_total {0, 0};
        set<Domain*> accounted_domains; // Domains whose energy is already counted
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                Occupancy state {domain->m_state};
                if (state == Occupancy::bound or state == Occupancy::misbound) {
                    if (find(accounted_domains.begin(), accounted_domains.end(),
                                domain) == accounted_domains.end()) {
                        Domain& bound_domain {*domain->m_bound_domain};
                        DH_DS_total.enthalpy += hybridization_enthalpy(*domain,
                                bound_domain);
                        DH_DS_total.entropy += hybridization_entropy(*domain,
                                bound_domain);
                        accounted_domains.insert(domain->m_bound_domain);
                    }
                }
            }
        }
        return DH_DS_total;
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

        // Check that energy has returned to 0 (eps is totally arbitrary)
        double eps {0.000001};
        if (m_energy < -eps or m_energy > eps) {
            cout << "Inconsistency in system energy\n";
            cout << m_energy << "\n";
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
        double delta_e {internal_check_domain_constraints(cd_i, pos, ore)};

        return delta_e;
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
        double delta_e {internal_unassign_domain(cd_i)};
        m_energy += delta_e;
        m_num_unassigned_domains++;

        return delta_e;
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
            size_t domain_size {m_sequences[c_i_ident][d_i].size()};
            if (domain_size == 16 or domain_size == 15) {
                domain = new SixteenDomain {c_i, c_i_ident, d_i, d_i_ident, chain_length};
            }
            else {
                throw NotImplemented {};
            }

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
            delta_e += hybridization_energy(cd_i, *cd_i.m_bound_domain);
        }
        else if (cd_i.m_state == Occupancy::bound) {
            delta_e += hybridization_energy(cd_i, *cd_i.m_bound_domain);
            delta_e += check_stacking(cd_i, *cd_i.m_bound_domain);
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
        double delta_e {internal_check_domain_constraints(cd_i, pos, ore)};
        if (not m_constraints_violated) {
            update_occupancies(cd_i, pos);
            m_energy += delta_e;
            m_num_unassigned_domains--;
        }
        return delta_e;
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

        // Update hybridization and stacking energy tables
        if (m_hybridization_energy_tables.count(temp) == 0) {
            get_energies();
            m_hybridization_energy_tables[temp] = m_hybridization_energies;
            m_hybridization_enthalpy_tables[temp] = m_hybridization_enthalpies;
            m_hybridization_entropy_tables[temp] = m_hybridization_entropies;
            m_stacking_energy_tables[temp] = m_stacking_energies;
        }
        else {
            unordered_map<pair<int, int>, double> m_hybrid_old {m_hybridization_energies};
            m_hybridization_energies = m_hybridization_energy_tables[temp];
            m_hybridization_enthalpies = m_hybridization_enthalpy_tables[temp];
            m_hybridization_entropies = m_hybridization_entropy_tables[temp];
            m_stacking_energies = m_stacking_energy_tables[temp];
        }

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

    // Protected methods

    double OrigamiSystem::bind_noncomplementary_domains(
            Domain& cd_i,
            Domain& cd_j) {
        m_constraints_violated = true;
        if (not check_domain_orientations_opposing(cd_i, cd_j)) {
            return 0;
        }
        m_constraints_violated = false;
        return OrigamiSystem::hybridization_energy(cd_i, cd_j);
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
                int scaffold_d_i {utility::index(m_identities[c_scaffold],
                        -staple_d_i)};
                scaffold_d_is.push_back(scaffold_d_i);
            }
            m_staple_ident_to_scaffold_ds.push_back(scaffold_d_is);
        }
    }

    void OrigamiSystem::get_energies() {
        // Get S, H, and G for all possible interactions and store
        if (m_energy_filebase.size() != 0) {
            if (not read_energies_from_file()) {
                calc_energies();
            }
        }
        else {
            calc_energies();
        }
    }

    void OrigamiSystem::calc_energies() {
        // Calculate S, H, and G for all possible interactions and store

        // Loop through all pairs of sequences
        for (size_t c_i {0}; c_i != m_sequences.size(); c_i++) {
            for (size_t c_j {0}; c_j != m_sequences.size(); c_j++) {
                size_t c_i_length {m_sequences[c_i].size()};
                size_t c_j_length {m_sequences[c_j].size()};
                for (size_t d_i {0}; d_i != c_i_length; d_i++) {
                    int d_i_ident {m_identities[c_i][d_i]};
                    for (size_t d_j {0}; d_j != c_j_length; d_j++) {
                        int d_j_ident {m_identities[c_j][d_j]};
                        string seq_i {m_sequences[c_i][d_i]};
                        string seq_j {m_sequences[c_j][d_j]};
                        pair<int, int> key {d_i_ident, d_j_ident};
                        calc_energy(seq_i, seq_j, key);
                    }
                }
            }
        }
        write_energies_to_file();
    }

    void OrigamiSystem::calc_energy(string seq_i, string seq_j,
            pair<int, int> key) {
        // Calculate S, H, and G for pair of sequences and store

        // Hybridization values
        vector<string> comp_seqs {nearestNeighbour::
                find_longest_contig_complement(seq_i, seq_j)};
        double H_hyb {0};
        double S_hyb {0};
        int N {0};
        
        // No interaction if no complementary sequence
        if (comp_seqs.size() == 0) {
            H_hyb = 0;
            S_hyb = 0;
        }

        // Take average value of H and S of all equal length comp seqs
        else {
            for (auto comp_seq: comp_seqs) {
                ThermoOfHybrid DH_DS {nearestNeighbour::
                        calc_unitless_hybridization_thermo(comp_seq,
                        m_temp, m_cation_M)};
                H_hyb += DH_DS.enthalpy;
                S_hyb += DH_DS.entropy;
                N++;
            }
            H_hyb /= N;
            S_hyb /= N;
        }
        m_hybridization_enthalpies[key] = H_hyb;
        m_hybridization_entropies[key] = S_hyb;
        m_hybridization_energies[key] = H_hyb - S_hyb;

        // Stacking energies
        double s_energy {0};
        s_energy += nearestNeighbour::calc_stacking_energy(seq_i, seq_j, m_temp,
                m_cation_M);
        m_stacking_energies[key] = s_energy;
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

    double OrigamiSystem::check_stacking(Domain& cd_new, Domain& cd_old) {
        // Check both sides
        double delta_e {0};
        for (int i: {-1, 0}) {
            Domain* cd_i {cd_old + i};
            Domain* cd_j {cd_old + (i + 1)};
            if (cd_i == nullptr or cd_j == nullptr) {
                continue;
            }
            bool cd_j_bound {cd_j->m_state == Occupancy::bound};
            if (not cd_j_bound) {
                continue;
            }
            Domain& cd_adjacent {*cd_j->m_bound_domain};
            if (cd_adjacent.m_c == cd_new.m_c) {
                continue;
            }
            VectorThree pos_i = {cd_i->m_pos};
            VectorThree pos_j = {cd_j->m_ore};
            VectorThree ndr {pos_j - pos_i};
            VectorThree ore_i = {cd_i->m_ore};
            if (ndr == ore_i) {
                continue;
            }
            delta_e += stacking_energy(cd_new, cd_adjacent);
        }
        return delta_e;
    }

    double OrigamiSystem::hybridization_energy(const Domain& cd_i,
            const Domain& cd_j) const {
        pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
        return m_hybridization_energies.at(key);
    }

    double OrigamiSystem::hybridization_enthalpy(const Domain& cd_i,
            const Domain& cd_j) const {
        pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
        return m_hybridization_enthalpies.at(key);
    }

    double OrigamiSystem::hybridization_entropy(const Domain& cd_i,
            const Domain& cd_j) const {
        pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
        return m_hybridization_entropies.at(key);
    }

    double OrigamiSystem::stacking_energy(const Domain& cd_i, const Domain& cd_j) const {
        pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
        return m_stacking_energies.at(key);
    }

    double OrigamiSystem::internal_unassign_domain(Domain& cd_i) {
        // Deletes positions, orientations, and removes/unassigns occupancies.
        Occupancy occupancy {cd_i.m_state};
        double delta_e {0};
        switch (occupancy) {
            case Occupancy::bound:
                m_num_fully_bound_domain_pairs -= 1;
                m_num_bound_domain_pairs -= 1;
                delta_e += check_stacking(cd_i, *cd_i.m_bound_domain);
                delta_e += unassign_bound_domain(cd_i);
                break;
            case Occupancy::misbound:
                m_num_bound_domain_pairs -= 1;
                if (cd_i.m_bound_domain->m_c == cd_i.m_c) {
                    m_num_self_bound_domain_pairs -= 1;
                }
                delta_e += unassign_bound_domain(cd_i);
                break;
            case Occupancy::unbound:
                unassign_unbound_domain(cd_i);
                break;
            case Occupancy::unassigned:

                // Allows for double unassignment
                m_num_unassigned_domains--;
                break;
        }
        return delta_e;
    }


    double OrigamiSystem::unassign_bound_domain(Domain& cd_i) {
        Domain& cd_j {*cd_i.m_bound_domain};
        double delta_e {-hybridization_energy(cd_i, cd_j)};

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

    bool OrigamiSystem::read_energies_from_file() {
        // Read energies from file, return false if not present
        string temp_string {"_" + std::to_string(static_cast<int>(m_temp))};

        // Hybridization free energies, enthalpies, entropies
        string henergy_filename {m_energy_filebase + temp_string + ".hene"};
        string hhenergy_filename {m_energy_filebase + temp_string + ".hhene"};
        string hsenergy_filename {m_energy_filebase + temp_string + ".hsene"};

        // Stacking energies
        string senergy_filename {m_energy_filebase + temp_string + ".sene"};

        std::ifstream henergy_file {henergy_filename};
        std::ifstream hhenergy_file {henergy_filename};
        std::ifstream hsenergy_file {henergy_filename};
        std::ifstream senergy_file {senergy_filename};

        bool files_present {henergy_file and hhenergy_file and hsenergy_file and
                senergy_file};
        if (files_present) {
            boost::archive::text_iarchive h_arch {henergy_file};
            h_arch >> m_hybridization_energies;
            boost::archive::text_iarchive hh_arch {hhenergy_file};
            hh_arch >> m_hybridization_enthalpies;
            boost::archive::text_iarchive hs_arch {hsenergy_file};
            hs_arch >> m_hybridization_entropies;
            boost::archive::text_iarchive s_arch {senergy_file};
            s_arch >> m_stacking_energies;
        }

        return files_present;
    }

    void OrigamiSystem::write_energies_to_file() {
        string temp_string {"_" + std::to_string(static_cast<int>(m_temp))};

        // Hybridization energies, enthalpies, entropies
        string henergy_filename {m_energy_filebase + temp_string + ".hene"};
        string hhenergy_filename {m_energy_filebase + temp_string + ".hhene"};
        string hsenergy_filename {m_energy_filebase + temp_string + ".hsene"};

        // Stacking energies
        string senergy_filename {m_energy_filebase + temp_string + ".sene"};

        std::ofstream henergy_file {henergy_filename};
        boost::archive::text_oarchive h_arch {henergy_file};
        h_arch << m_hybridization_energies;
        std::ofstream hhenergy_file {hhenergy_filename};
        boost::archive::text_oarchive hh_arch {hhenergy_file};
        hh_arch << m_hybridization_enthalpies;
        std::ofstream hsenergy_file {hsenergy_filename};
        boost::archive::text_oarchive hs_arch {hsenergy_file};
        hs_arch << m_hybridization_entropies;
        std::ofstream senergy_file {senergy_filename};
        boost::archive::text_oarchive s_arch {senergy_file};
        s_arch << m_stacking_energies;
    }

    void OrigamiSystem::update_energy() {
        // Recalculate total system energy with new energy tables
        // WARNING: Does not account for stacking energy
        m_energy = 0;
        set<Domain*> accounted_domains; // Domains whose energy is already counted
        for (auto chain: m_domains) {
            for (auto domain: chain) {
                Occupancy state {domain->m_state};
                if (state == Occupancy::bound or state == Occupancy::misbound) {
                    if (find(accounted_domains.begin(), accounted_domains.end(),
                                domain) == accounted_domains.end()) {
                        Domain& bound_domain {*domain->m_bound_domain};
                        m_energy += hybridization_energy(*domain, bound_domain);
                        accounted_domains.insert(domain->m_bound_domain);
                    }
                }
            }
        }
    }

    double OrigamiSystem::internal_check_domain_constraints(
            Domain& cd_i,
            VectorThree pos,
            VectorThree ore) {
        // Updates positions and orientations and returns without reverting if no 
        // constraint violation. But states left unassigned.
        Occupancy occupancy {position_occupancy(pos)};
        double delta_e {0};

        switch (occupancy) {
            case Occupancy::bound:
                m_constraints_violated = true;
                break;
            case Occupancy::misbound:
                m_constraints_violated = true;
                break;
            case Occupancy::unbound:
                update_domain(cd_i, pos, ore);
                update_occupancies(cd_i, pos);
                delta_e += bind_domain(cd_i);
                internal_unassign_domain(cd_i);
                break;
            case Occupancy::unassigned:
                update_domain(cd_i, pos, ore);
        }
        return delta_e;
    }

    double OrigamiSystem::bind_domain(Domain& cd_i) {
        // Check constraints for an unbound to (mis)bound transition
        Domain& cd_j {*cd_i.m_bound_domain};
        double delta_e {0};
        bool comp {check_domains_complementary(cd_i, cd_j)};
        if (comp) {
            delta_e += bind_complementary_domains(cd_i, cd_j);
        }
        else {
            delta_e += bind_noncomplementary_domains(cd_i, cd_j);
        }
        return delta_e;
    }

    double OrigamiSystem::bind_complementary_domains(Domain& cd_i, Domain& cd_j) {
        // cd_i is new
        m_constraints_violated = true;
        if (not check_domain_orientations_opposing(cd_i, cd_j)) {
            return 0;
        }
        if (not check_domain_pair_constraints(cd_i)) {
            return 0;
        }
        if (not check_domain_pair_constraints(cd_j)) {
            return 0;
        }

        // Missed one linear helix check per chain
        if (not check_linear_helix_rear(cd_i)) {
            return 0;
        }
        if (not check_linear_helix_rear(cd_j)) {
            return 0;
        }

        // Missed two contiguous junction checks per chain
        if (not check_junction_front(cd_i)) {
            return 0;
        }
        if (not check_junction_front(cd_j)) {
            return 0;
        }
        if (not check_junction_rear(cd_i)) {
            return 0;
        }
        if (not check_junction_rear(cd_j)) {
            return 0;
        }

        // Collect energies
        m_constraints_violated = false;
        double delta_e {0};
        delta_e += hybridization_energy(cd_i, cd_j);
        delta_e += check_stacking(cd_i, cd_j);

        return delta_e;
    }

    bool OrigamiSystem::check_domain_pair_constraints(Domain& cd_i) {
        bool pair_constraints_obeyed {true};
        // Check both pairs
        for (int i: {-1, 0}) {
            Domain* cd_1 {cd_i + i};
            Domain* cd_2 {cd_i + (i + 1)};
            if (cd_1 == nullptr or cd_2 == nullptr) {
                continue;
            }
            if (cd_1->m_state == Occupancy::bound and cd_2->m_state == Occupancy::bound) {
                if (not check_helical_constraints(*cd_1, *cd_2)) {
                    pair_constraints_obeyed = false;
                    break;
                }
            }
            else {
            }
        }
        return pair_constraints_obeyed;
    }

    bool OrigamiSystem::check_helical_constraints(Domain& cd_1, Domain& cd_2) {
        // Given that d1 and d2 exist and are bound
        bool helical_constraints_obeyed {true};

        // Calculate next domain vector
        VectorThree ndr {cd_2.m_pos - cd_1.m_pos};

        Domain& cd_bound_1 {*cd_1.m_bound_domain};
        Domain& cd_bound_2 {*cd_2.m_bound_domain};
        bool bound_same_chain {cd_bound_1.m_c == cd_bound_2.m_c};

        // New helix case
        if (cd_1.m_ore == ndr) {

            // Check doubly contiguous constraint
            if (bound_same_chain) {
                if (cd_bound_1.m_d == cd_bound_2.m_d - 1) {
                    helical_constraints_obeyed = false;
                }
                else {
                    if (not check_doubly_contiguous_junction(cd_1, cd_2)) {
                        helical_constraints_obeyed = false;
                    }
                }
            }
        }

        // Non-physical case
        else if (cd_1.m_ore == -ndr) {
            helical_constraints_obeyed = false;
        }

        // Same helix case
        else {

            // Check double contiguous constraint
            if (bound_same_chain) {
                if (cd_bound_1.m_d == cd_bound_2.m_d + 1) {
                    helical_constraints_obeyed = false;
                }
            }
            if (cd_1.check_twist_constraint(ndr, cd_2)) {
                if (not check_linear_helix(ndr, cd_2)) {
                    helical_constraints_obeyed = false;
                }
            }
            else {
                helical_constraints_obeyed = false;
            }
        }
        return helical_constraints_obeyed;
    }

    bool OrigamiSystem::check_domain_orientations_opposing(Domain& cd_i, Domain& cd_j) {
        bool domain_orientations_opposing {true};
        if (cd_i.m_ore != -cd_j.m_ore) {
            domain_orientations_opposing = false;
            return domain_orientations_opposing;
        }
        return domain_orientations_opposing;
    }

    bool OrigamiSystem::check_linear_helix(VectorThree ndr_1,Domain& cd_2) {
        bool linear_if_helix {true};
        Domain* cd_3 {cd_2 + 1};
        if (cd_3 == nullptr) {
            return linear_if_helix;
        }

        // Third domain not bound or unasssigned
        if (cd_3->m_state != Occupancy::bound) {
            return linear_if_helix;
        }

        // Third domain part of new helix
        VectorThree ndr_2 {cd_3->m_pos - cd_2.m_pos};
        if (ndr_2 == cd_2.m_ore) {
            return linear_if_helix;
        }

        // Third domain linear
        if (ndr_1 == ndr_2) {
            return linear_if_helix;
        }

        linear_if_helix = false;
        return linear_if_helix;
    }

    bool OrigamiSystem::check_linear_helix_rear(Domain& cd_3) {
        // Check linear helix constraints given rear domain exists and bound
        bool linear_if_helix {true};
        Domain* cd_2 {cd_3 + -1};
        Domain* cd_1 {cd_3 + -2};
        if (cd_1 == nullptr or cd_2 == nullptr) {
            return linear_if_helix;
        }
        
        // Domains 1 and 2 not bound
        bool cd_1_bound {cd_1->m_state == Occupancy::bound};
        bool cd_2_bound {cd_2->m_state == Occupancy::bound};
        if (not (cd_1_bound and cd_2_bound)) {
            return linear_if_helix;
        }
        VectorThree ndr_1 {cd_2->m_pos - cd_1->m_pos};
        VectorThree ndr_2 {cd_3.m_pos - cd_2->m_pos};
        VectorThree ore_1 {cd_1->m_ore};
        VectorThree ore_2 {cd_2->m_ore};
        // Domains 1 and 2 are junction
        if (ndr_1 == ore_1) {
            return linear_if_helix;
        }

        // Domains 2 and 3 are junction
        if (ndr_2 == ore_2) {
            return linear_if_helix;
        }

        // Domains are linear
        if (ndr_1 == ndr_2) {
            return linear_if_helix;
        }

        linear_if_helix = false;
        return linear_if_helix;
    }

    bool OrigamiSystem::check_doubly_contiguous_junction(Domain& cd_2, Domain& cd_3) {
        // Already know d_2 and d_3 are doubly contiguous junction
        bool junction_constraints_obeyed {true};
        Domain* cd_1 {cd_2 + -1};
        Domain* cd_4 {cd_3 + 1};
        if (cd_1 == nullptr or cd_4 == nullptr) {
            return junction_constraints_obeyed;
        }

        bool cd_1_bound {cd_1->m_state == Occupancy::bound};
        bool cd_4_bound {cd_4->m_state == Occupancy::bound};
        if (not (cd_1_bound and cd_4_bound)) {
            return junction_constraints_obeyed;
        }
        junction_constraints_obeyed = check_doubly_contiguous_junction(*cd_1, cd_2,
                cd_3, *cd_4);
        return junction_constraints_obeyed;
    }

    bool OrigamiSystem::doubly_contiguous_junction(Domain& cd_1, Domain& cd_2) {
        // Given that d1 and d2 exist and are bound and contiguous

        // Check doubly contiguous
        Domain& cd_bound_1 {*cd_1.m_bound_domain};
        Domain& cd_bound_2 {*cd_2.m_bound_domain};
        if (cd_bound_1.m_c != cd_bound_2.m_c) {
            return false;
        }

        if (cd_bound_2.m_d != cd_bound_1.m_d - 1) {
            return false;
        }

        // Check junction
        VectorThree ndr {cd_2.m_pos - cd_1.m_pos};
        VectorThree ore_1 {cd_1.m_ore};
        if (ndr != ore_1) {

            // To be clear, this is only used to check for changes in junction
            // constraints; if it's been set before the local constraints were okay
            return false;
        }

        return true;
    }

    bool OrigamiSystem::check_doubly_contiguous_junction(
            Domain& cd_1,
            Domain& cd_2,
            Domain& cd_3,
            Domain& cd_4) {
        // Calling functions already check for existance and boundeness of d1-4.

        bool junction_constraints_obeyed;
        VectorThree ndr_1 {cd_2.m_pos - cd_1.m_pos};
        VectorThree ndr_3 {cd_4.m_pos - cd_3.m_pos};
        if (ndr_1 == -ndr_3) {
        // if (ndr_1 != ndr_3) {
            junction_constraints_obeyed = true;
        }
        else {
            junction_constraints_obeyed = false;
        }

        return junction_constraints_obeyed;
    }

    bool OrigamiSystem::check_junction_front(Domain& cd_1) {
        bool junction_constraints_obeyed {true};
        Domain* cd_2 {cd_1 + 1};
        Domain* cd_3 {cd_1 + 2};
        Domain* cd_4 {cd_1 + 3};
        if (cd_2 == nullptr or cd_3 == nullptr or cd_4 == nullptr) {
            return junction_constraints_obeyed;
        }

        // Check that all domains in bound state
        bool cd_2_bound {cd_2->m_state == Occupancy::bound};
        bool cd_3_bound {cd_3->m_state == Occupancy::bound};
        bool cd_4_bound {cd_4->m_state == Occupancy::bound};
        if (not (cd_2_bound and cd_3_bound and cd_4_bound)) {
            return junction_constraints_obeyed;
        }

        if (not doubly_contiguous_junction(*cd_2, *cd_3)) {
            return junction_constraints_obeyed;
        }

        return check_doubly_contiguous_junction(cd_1, *cd_2, *cd_3, *cd_4);
    }

    bool OrigamiSystem::check_junction_rear(Domain& cd_4) {
        // Check junction given d4 exists
        bool junction_constraints_obeyed {true};
        Domain* cd_3 {cd_4 + -1};
        Domain* cd_2 {cd_4 + -2};
        Domain* cd_1 {cd_4 + -3};
        if (cd_1 == nullptr or cd_2 == nullptr or cd_3 == nullptr) {
            return junction_constraints_obeyed;
        }

        // Check that all domains in bound state
        bool cd_1_bound {cd_1->m_state == Occupancy::bound};
        bool cd_2_bound {cd_2->m_state == Occupancy::bound};
        bool cd_3_bound {cd_3->m_state == Occupancy::bound};
        if (not (cd_1_bound and cd_2_bound and cd_3_bound)) {
            return junction_constraints_obeyed;
        }

        if (not doubly_contiguous_junction(*cd_2, *cd_3)) {
            return junction_constraints_obeyed;
        }

        return check_doubly_contiguous_junction(*cd_1, *cd_2, *cd_3, cd_4);
    }

    double OrigamiSystemWithoutMisbinding::bind_noncomplementary_domains(
            Domain&, Domain&) {
        m_constraints_violated = true;
        return 0;
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
            if (check_domains_complementary(*unbound_domain_at(pos), cd_i)) {
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

    double molarity_to_lattice_volume(double molarity, double lattice_site_volume) {
        // Given a molarity, calculate the volume that cancels the fugacity.
        // Volume is in units of number of lattice sites.

        // Number of lattice sites per L (1 L * (1000) cm^3 / L * m^3 / (10^2)^3 cm^3)
        double sites_per_litre {1e-3 / lattice_site_volume};

        // u = KB*T*ln(p), where p is the number of particles per lattice site
        // z = exp(1/(KB*T)*u) = exp(ln(p)) = p
        // V * p = 1, V = 1 / p
        // So just convert molarity to number of particles per lattice site
        double V {1 / (molarity * utility::NA / sites_per_litre)};
        return V;
    }

    double molarity_to_chempot(double molarity, double temp,
            double lattice_site_volume) {
        double sites_per_litre {1e-3 / lattice_site_volume};
        double chempot {temp * log(molarity * utility::NA / sites_per_litre)};
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
                params.m_temp_for_staple_u, params.m_lattice_site_volume)};
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
