// origami_system.cpp

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/utility.hpp>

#include "origami_system.h"
#include "utility.h"
#include "nearest_neighbour.h"

using std::abs;
using std::max_element;
using std::cout;

using namespace NearestNeighbour;
using namespace Utility;
using namespace Origami;

// Public methods

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
        double temp,
        double lattice_site_volume,
        double cation_M,
        double staple_M,
        bool cyclic,
        string energy_filebase) :

        m_identities {identities},
        m_sequences {sequences},
        m_temp {temp},
        m_cation_M {cation_M},
        m_staple_M {staple_M},
    	m_volume {molarity_to_lattice_volume(m_staple_M, lattice_site_volume)},
        m_cyclic {cyclic},
        m_energy_filebase {energy_filebase} {

    initialize_complementary_associations();
    initialize_config(chains);
    update_temp(m_temp);
}

OrigamiSystem::~OrigamiSystem() {
    for (auto chain: m_domains) {
        for (auto domain: chain) {
            delete domain;
        }
    }
}

vector<Domain*> OrigamiSystem::get_chain(int c_i) {
    int c_i_index {index(m_chain_indices, c_i)};
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

double OrigamiSystem::unassign_domain(Domain& cd_i) {
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
            break;
    }
    m_energy += delta_e;
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

void OrigamiSystem::centre() {
    // Translate the system such that the first scaffold domain is on the origin

    // The maps that go from a position to something need to be reset
    unordered_map<VectorThree, Domain*> pos_to_unbound_d {};
    unordered_map<VectorThree, Occupancy> position_occupancies {};
    VectorThree refpos {m_domains[0][0]->m_pos};
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

void OrigamiSystem::update_temp(double temp) {
    m_temp = temp;
    if (m_hybridization_energy_tables.count(temp) == 0) {
        if (m_energy_filebase.size() != 0) {
            read_energies_from_file(temp);
        }
        else {
            calculate_energies();
        }
        m_hybridization_energy_tables[temp] = m_hybridization_energies;
        m_stacking_energy_tables[temp] = m_stacking_energies;
    }
    else {
        unordered_map<pair<int, int>, double> m_hybrid_old {m_hybridization_energies};
        m_hybridization_energies = m_hybridization_energy_tables[temp];
        m_stacking_energies = m_stacking_energy_tables[temp];
    }
    update_energy();
}

void OrigamiSystem::read_energies_from_file(double temp) {
    string temp_string {"_" + std::to_string(static_cast<int>(temp))};
    string henergy_filename {m_energy_filebase + temp_string + ".hene"};
    string senergy_filename {m_energy_filebase + temp_string + ".sene"};
    std::ifstream henergy_file {henergy_filename};
    std::ifstream senergy_file {senergy_filename};
    if (henergy_file and senergy_file) {
        boost::archive::text_iarchive h_arch {henergy_file};
        h_arch >> m_hybridization_energies;
        boost::archive::text_iarchive s_arch {senergy_file};
        s_arch >> m_stacking_energies;
    }
    else {
        calculate_energies();
        std::ofstream henergy_file {henergy_filename};
        boost::archive::text_oarchive h_arch {henergy_file};
        h_arch << m_hybridization_energies;
        std::ofstream senergy_file {senergy_filename};
        boost::archive::text_oarchive s_arch {senergy_file};
        s_arch << m_stacking_energies;
    }
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

double OrigamiSystem::set_domain_config(
        Domain& cd_i,
        VectorThree pos,
        VectorThree ore) {
    if (cd_i.m_state != Occupancy::unassigned) {
        cout << "e\n";
        throw OrigamiMisuse {};
    }

    // Check constraints and update if obeyed, otherwise throw
    double delta_e {check_domain_constraints(cd_i, pos, ore)};
    if (not m_constraints_violated) {
        update_occupancies(cd_i, pos);
        m_energy += delta_e;
    }
    return delta_e;
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
    return delta_e;
}

double OrigamiSystem::check_domain_constraints(
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

        m_num_domains += 1;
    }

    return c_i;
}

void OrigamiSystem::delete_chain(int c_i) {
    // Delete chain c_i_
    
    int c_i_index {index(m_chain_indices, c_i)};
    int c_i_ident {m_chain_identities[c_i_index]};

    // Index in m_identity_to_index of given index and type
    int j {index(m_identity_to_index[c_i_ident], c_i)};
    m_identity_to_index[c_i_ident].erase(m_identity_to_index[c_i_ident].begin()
            + j);
    m_chain_indices.erase(m_chain_indices.begin() + c_i_index);
    m_chain_identities.erase(m_chain_identities.begin() + c_i_index);
    m_num_domains -= m_domains[c_i_index].size();
    for (auto domain: m_domains[c_i_index]) {
        delete domain;
    }
    m_domains.erase(m_domains.begin() + c_i_index);
}

// Private methods

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
            int scaffold_d_i {index(m_identities[c_scaffold], -staple_d_i)};
            scaffold_d_is.push_back(scaffold_d_i);
        }
        m_staple_ident_to_scaffold_ds.push_back(scaffold_d_is);
    }
}

void OrigamiSystem::initialize_config(Chains chains) {
    // Must still run check_all_constraints

    // Create domain objects
    for (size_t i {0}; i != chains.size(); i++) {
        Chain chain {chains[i]};
        int c_i {chain.index};
        int c_i_ident {chain.identity};
        add_chain(c_i_ident, c_i);

        // Make scaffold chain domains modular if cyclic
        if (m_cyclic and c_i == c_scaffold) {
            Domain* first_domain {m_domains[c_i][0]};
            Domain* last_domain {m_domains[c_i].back()};
            last_domain->m_forward_domain = first_domain;
            first_domain->m_backward_domain = last_domain;
        }
    }

    // Set configuration of domains
    for (size_t i {0}; i != chains.size(); i++) {
        Chain chain {chains[i]};
        int c_i {chain.index};
        int domains_c_i {m_chain_indices[c_i]};
        int num_domains {static_cast<int>(m_domains[domains_c_i].size())};
        for (int d_i {0}; d_i != num_domains; d_i++) {
            Domain* domain {m_domains[c_i][d_i]};
            VectorThree pos = chain.positions[d_i];
            VectorThree ore = chain.orientations[d_i];
            domain->m_pos = pos;
            domain->m_ore = ore;
        }
    }

    m_current_c_i = *max_element(m_chain_identities.begin(),
            m_chain_identities.end());
}

void OrigamiSystem::calculate_energies() {
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

                    // Hybridization energies
                    vector<string> comp_seqs {find_longest_contig_complement(
                            seq_i, seq_j)};
                    double energy {0};
                    int N {0};
                    if (comp_seqs.size() == 0) {
                        energy = 0;
                    }
                    else {
                        for (auto comp_seq: comp_seqs) {
                            energy += calc_hybridization_energy(
                                   comp_seq, m_temp, m_cation_M);
                            N++;
                        }
                        energy /= N;
                    }
                    m_hybridization_energies[key] = energy;

                    // Stacking energies
                    double s_energy {0};
                    s_energy += calc_stacking_energy(seq_i, seq_j, m_temp,
                            m_cation_M);
                    m_stacking_energies[key] = s_energy;
                }
            }
        }
    }
}

double OrigamiSystem::hybridization_energy(const Domain& cd_i,
        const Domain& cd_j) const {
    pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
    return m_hybridization_energies.at(key);
}

double OrigamiSystem::stacking_energy(const Domain& cd_i, const Domain& cd_j) const {
    pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
    return m_stacking_energies.at(key);
}

void OrigamiSystem::internal_unassign_domain(Domain& cd_i) {
    m_energy -= unassign_domain(cd_i);
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
            cout << "a\n";
            throw OrigamiMisuse {};
    }
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

double OrigamiSystem::bind_noncomplementary_domains(
        Domain& cd_i,
        Domain& cd_j) {
    return OrigamiSystem::hybridization_energy(cd_i, cd_j);
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

void OrigamiSystem::update_energy() {

    // Unassign everything (and check nothing was already unassigned)
    for (auto chain: m_domains) {
        for (auto domain: chain) {
            unassign_domain(*domain);
        }
    }
    m_energy = 0;

    // Reset configuration
    for (auto chain: m_domains) {
        for (auto domain: chain) {
            set_domain_config(*domain, domain->m_pos, domain->m_ore);
        }
    }
}

void OrigamiSystem::check_all_constraints() {

    // Unassign everything (and check nothing was already unassigned)
    for (auto chain: m_domains) {
        for (auto domain: chain) {
            if (domain->m_state == Occupancy::unassigned) {
                cout << "g\n";
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
        throw OrigamiMisuse {};
    }
    else {

        // Prevent rounding errors from building up
        m_energy = 0;
    }

    // Reset configuration
    for (auto chain: m_domains) {
        for (auto domain: chain) {
            set_domain_config(*domain, domain->m_pos, domain->m_ore);
            if (m_constraints_violated) {
                cout << "b\n";
                set_domain_config(*domain, domain->m_pos, domain->m_ore);
                throw OrigamiMisuse {};
            }
        }
    }
    check_distance_constraints();
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
                    cout << "c\n";
                    throw OrigamiMisuse {};
                }
            }
        }
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
        junction_constraints_obeyed = true;
    }
    else {
        junction_constraints_obeyed = false;
    }

    return junction_constraints_obeyed;
}

bool OrigamiSystem::check_domain_orientations_opposing(Domain& cd_i, Domain& cd_j) {
    bool domain_orientations_opposing {true};
    if (cd_i.m_ore != -cd_j.m_ore) {
        domain_orientations_opposing = false;
        return domain_orientations_opposing;
    }
    return domain_orientations_opposing;
}

double OrigamiSystemWithoutMisbinding::bind_noncomplementary_domains(
        Domain&, Domain&) {
    m_constraints_violated = true;
    return 0;
}

double Origami::molarity_to_lattice_volume(double molarity, double lattice_site_volume) {
    // Given a molarity, calculate the volume that cancels the fugacity.
    // Volume is in units of number of lattice sites.

    // Number of lattice sites per L (1 L * (1000) cm^3 / L * m^3 / (10^2)^3 cm^3)
    double sites_per_litre {1e-3 / lattice_site_volume};

    // u = KB*T*ln(p), where p is the number of particles per lattice site
    // g = exp(1/(-KB*T)*u) = exp(ln(p)) = p
    // V * p = 1, V = 1 / p
    // So just convert molarity to number of particles per lattice site
    double V {1 / (molarity * Utility::NA / sites_per_litre)};
    return V;
}
