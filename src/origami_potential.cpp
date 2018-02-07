// origami_potential.cpp

#include "origami_potential.h"

#include <iostream>
#include <fstream>
#include <utility>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/utility.hpp>

namespace potential {

    using utility::Occupancy;

    bool check_domain_orientations_opposing(Domain& cd_i, Domain& cd_j) {
        bool domain_orientations_opposing {true};
        if (cd_i.m_ore != -cd_j.m_ore) {
            domain_orientations_opposing = false;
            return domain_orientations_opposing;
        }
        return domain_orientations_opposing;
    }

    bool check_domains_exist_and_bound(vector<Domain*> cdv) {
        bool exists_and_bound {true};
        for (auto cd: cdv) {
            if (cd == nullptr) {
                exists_and_bound = false;
                break;
            }
            bool cd_bound {cd->m_state == Occupancy::bound};
            if (not cd_bound) {
                exists_and_bound = false;
                break;
            }
        }

        return exists_and_bound;
    }

    BindingPotential::BindingPotential(OrigamiPotential& pot):
            m_pot {pot} {
    }

    pair<double, int> BindingPotential::check_stacking(Domain& cd_i, Domain& cd_j) {

        // Loop through relevant pairs
        double delta_e {0};
        int stacked_pairs {0};
        for (auto cd: {&cd_i, &cd_j}) {
            for (int i: {-1, 0}) {
                Domain* cd_1 {*cd + i};
                Domain* cd_2 {*cd + (i + 1)};
                if (not check_domains_exist_and_bound({cd_1, cd_2})) {
                    continue;
                }
                Domain& cd_bound_1 {*cd_1->m_bound_domain};
                Domain& cd_bound_2 {*cd_2->m_bound_domain};
                bool bound_same_chain {cd_bound_1.m_c == cd_bound_2.m_c};
                double delta_delta_e;
                int delta_stacked_pairs;
                if (bound_same_chain and
                        (cd_bound_1.m_d == cd_bound_2.m_d - 1)) {

                    // I have to divide this by to so that I don't overcount
                    std::tie(delta_delta_e, delta_stacked_pairs) =
                        check_pair_stacking(cd_1, cd_2);
                    delta_e += delta_delta_e/2;
                    //DEBUG
                    //stacked_pairs += 1;
                }
                else {
                    std::tie(delta_delta_e, delta_stacked_pairs) =
                        check_pair_stacking(cd_1, cd_2);
                    delta_e += delta_delta_e;
                    stacked_pairs += delta_stacked_pairs;
                }
            }
        }

        return {delta_e, stacked_pairs};
    }

    pair<double, int> RestrictiveBindingPotential::bind_domains(
            Domain& cd_i, Domain& cd_j) {

        // cd_i is new
        m_constraints_violated = true;
        if (not check_domain_orientations_opposing(cd_i, cd_j)) {
            return {0, 0};
        }
        if (not check_domain_pair_constraints(cd_i)) {
            return {0, 0};
        }
        if (not check_domain_pair_constraints(cd_j)) {
            return {0, 0};
        }

        // Missed one linear helix check per chain
        if (not check_linear_helix_rear(cd_i)) {
            return {0, 0};
        }
        if (not check_linear_helix_rear(cd_j)) {
            return {0, 0};
        }

        // Missed two contiguous junction checks per chain
        if (not check_junction_front(cd_i)) {
            return {0, 0};
        }
        if (not check_junction_front(cd_j)) {
            return {0, 0};
        }
        if (not check_junction_rear(cd_i)) {
            return {0, 0};
        }
        if (not check_junction_rear(cd_j)) {
            return {0, 0};
        }

        // Collect energies
        m_constraints_violated = false;
        double delta_e {0};
        int stacked_pairs {0};
        std::tie(delta_e, stacked_pairs) = check_stacking(cd_i, cd_j);
        delta_e += m_pot.hybridization_energy(cd_i, cd_j);

        return {delta_e, stacked_pairs};
    }

    bool RestrictiveBindingPotential::check_domain_pair_constraints(
            Domain& cd_i) {

        // Check both pairs
        bool pair_constraints_obeyed {true};
        for (int i: {-1, 0}) {
            Domain* cd_1 {cd_i + i};
            Domain* cd_2 {cd_i + (i + 1)};
            if (not check_domains_exist_and_bound({cd_1, cd_2})) {
                continue;
            }
            if (not check_helical_constraints(*cd_1, *cd_2)) {
                pair_constraints_obeyed = false;
                break;
            }
        }

        return pair_constraints_obeyed;
    }

    bool RestrictiveBindingPotential::check_helical_constraints(Domain& cd_1,
            Domain& cd_2) {

        // Given that d1 and d2 exist and are bound
        bool helical_constraints_obeyed {true};

        // Calculate next domain vector
        VectorThree ndr {cd_2.m_pos - cd_1.m_pos};

        Domain& cd_bound_1 {*cd_1.m_bound_domain};
        Domain& cd_bound_2 {*cd_2.m_bound_domain};
        bool bound_same_chain {cd_bound_1.m_c == cd_bound_2.m_c};

        // New helix case
        if (cd_1.m_ore == ndr) {

            // Check parallel helix constraint
            if (not (cd_1.m_ore == cd_2.m_ore)) {
                helical_constraints_obeyed = false;
            }

            // Check doubly contiguous constraint
            if (helical_constraints_obeyed and bound_same_chain) {
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

    bool RestrictiveBindingPotential::check_linear_helix(VectorThree ndr_1,
            Domain& cd_2) {

        bool linear_if_helix {true};
        Domain* cd_3 {cd_2 + 1};
        if (not check_domains_exist_and_bound({cd_3})) {
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

    bool RestrictiveBindingPotential::check_linear_helix_rear(Domain& cd_3) {
        // Check linear helix constraints given rear domain exists and bound
        bool linear_if_helix {true};
        Domain* cd_2 {cd_3 + -1};
        Domain* cd_1 {cd_3 + -2};
        if (not check_domains_exist_and_bound({cd_1, cd_2})) {
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

    bool RestrictiveBindingPotential::check_doubly_contiguous_junction(
            Domain& cd_2, Domain& cd_3) {

        // Already know d_2 and d_3 are doubly contiguous junction
        bool junction_constraints_obeyed {true};
        Domain* cd_1 {cd_2 + -1};
        Domain* cd_4 {cd_3 + 1};
        if (not check_domains_exist_and_bound({cd_1, cd_4})) {
            return junction_constraints_obeyed;
        }
        junction_constraints_obeyed = check_doubly_contiguous_junction(*cd_1, cd_2,
                cd_3, *cd_4);
        return junction_constraints_obeyed;
    }

    bool RestrictiveBindingPotential::doubly_contiguous_junction(Domain& cd_1,
            Domain& cd_2) {

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

    bool RestrictiveBindingPotential::check_doubly_contiguous_junction(
            Domain& cd_1,
            Domain& cd_2,
            Domain& cd_3,
            Domain& cd_4) {

        bool junction_constraints_obeyed;
        VectorThree ndr_1 {cd_2.m_pos - cd_1.m_pos};
        VectorThree ndr_3 {cd_4.m_pos - cd_3.m_pos};

        if (ndr_1 == -ndr_3 or ndr_1 == cd_1.m_ore or ndr_3 == cd_3.m_ore) {
            junction_constraints_obeyed = true;
        }
        else {
            junction_constraints_obeyed = false;
        }

        return junction_constraints_obeyed;
    }

    bool RestrictiveBindingPotential::check_junction_front(Domain& cd_1) {
        bool junction_constraints_obeyed {true};
        Domain* cd_2 {cd_1 + 1};
        Domain* cd_3 {cd_1 + 2};
        Domain* cd_4 {cd_1 + 3};
        if (not check_domains_exist_and_bound({cd_2, cd_3, cd_4})) {
            return junction_constraints_obeyed;
        }

        if (not doubly_contiguous_junction(*cd_2, *cd_3)) {
            return junction_constraints_obeyed;
        }

        return check_doubly_contiguous_junction(cd_1, *cd_2, *cd_3, *cd_4);
    }

    bool RestrictiveBindingPotential::check_junction_rear(Domain& cd_4) {
        // Check junction given d4 exists
        bool junction_constraints_obeyed {true};
        Domain* cd_3 {cd_4 + -1};
        Domain* cd_2 {cd_4 + -2};
        Domain* cd_1 {cd_4 + -3};
        if (not check_domains_exist_and_bound({cd_1, cd_2, cd_3})) {
            return junction_constraints_obeyed;
        }

        if (not doubly_contiguous_junction(*cd_2, *cd_3)) {
            return junction_constraints_obeyed;
        }

        return check_doubly_contiguous_junction(*cd_1, *cd_2, *cd_3, cd_4);
    }

    pair<double,int> RestrictiveBindingPotential::check_pair_stacking(Domain* cd_1,
            Domain* cd_2) {

        double delta_e {0};
        int stacked {0};
        VectorThree ndr {cd_2->m_pos - cd_1->m_pos};
        if (ndr == cd_1->m_ore) {
            ;
        }
        else {
            delta_e += m_pot.stacking_energy(*cd_1, *cd_2);
            stacked = 1;
        }

        return {delta_e, stacked};
    }

    pair<double, int> FlexibleBindingPotential::bind_domains(
            Domain& cd_i,
            Domain& cd_j) {

        m_constraints_violated = false;
        if (not check_domain_orientations_opposing(cd_i, cd_j)) {
            m_constraints_violated = true;
            return {0, 0};
        }

        // Loop through relevant pairs
        double delta_e {0};
        int stacked_pairs {0};
        pair<double, int> delta_ene_stack;
        for (auto cd: {&cd_i, &cd_j}) {
            for (int i: {-1, 0}) {
                Domain* cd_1 {*cd + i};
                Domain* cd_2 {*cd + (i + 1)};
                if (not check_domains_exist_and_bound({cd_1, cd_2})) {
                    continue;
                }
                delta_ene_stack = check_constraints(cd_1, cd_2, i);
                delta_e += delta_ene_stack.first;
                stacked_pairs += delta_ene_stack.second;
                if (m_constraints_violated) {
                    return {0, 0};
                }
            }
        }
        delta_e += m_pot.hybridization_energy(cd_i, cd_j);

        return {delta_e, stacked_pairs};
    }

    pair<double, int> FlexibleBindingPotential::check_constraints(
            Domain* cd_1, Domain* cd_2, int i) {

        // This assumes that the staples are at most 2 domains long
        // Othewise I would have to check the edge pair regardless of
        // whether the core was doubly contiguous
        double delta_e {0};
        int stacked_pairs {0};
        Domain& cd_bound_1 {*cd_1->m_bound_domain};
        Domain& cd_bound_2 {*cd_2->m_bound_domain};
        bool bound_same_chain {cd_bound_1.m_c == cd_bound_2.m_c};
        VectorThree ndr {cd_2->m_pos - cd_1->m_pos};
        pair<double, int> delta_ene_stack;
        if (bound_same_chain) {

            // Same helix case
            if (cd_bound_1.m_d == cd_bound_2.m_d - 1) {
                if (ndr == cd_1->m_ore or ndr == -cd_1->m_ore) {
                    m_constraints_violated = true;
                }
                else if (cd_1->check_twist_constraint(ndr, *cd_2)) {

                    // I have to divide this by to so that I don't overcount
                    delta_e += m_pot.stacking_energy(*cd_1, *cd_2) / 2;
                    //DEBUG
                    //stacked_pairs += 1;
                }
                else {
                    m_constraints_violated = true;
                }
            }

            // Crossover case
            else if (cd_bound_1.m_d == cd_bound_2.m_d + 1) {
                if (cd_1->m_ore == ndr) {
                    Domain* cd_j1 {*cd_1 + -1};
                    Domain* cd_j4 {*cd_1 + -1};
                    if (check_domains_exist_and_bound({cd_j1, cd_j4})) {
                        if (not check_doubly_contig_junction(cd_j1, cd_1, cd_2, cd_j4)) {
                            m_constraints_violated = true;
                        }
                    }
                }
                else {
                    m_constraints_violated = true;
                }
            }

            // ugly
            else {
                if (not check_edge_pair_junction(cd_1, cd_2, i)) {
                    m_constraints_violated = true;
                }
                delta_ene_stack = check_pair_stacking(cd_1, cd_2);
                delta_e += delta_ene_stack.first;
                stacked_pairs += delta_ene_stack.second;
            }
        }
        else {
            if (not check_edge_pair_junction(cd_1, cd_2, i)) {
                m_constraints_violated = true;
            }
            delta_ene_stack = check_pair_stacking(cd_1, cd_2);
            delta_e += delta_ene_stack.first;
            stacked_pairs += delta_ene_stack.second;
        }

        return {delta_e, stacked_pairs};
    }

    bool FlexibleBindingPotential::check_doubly_contig_junction(
            Domain* cd_1, Domain* cd_2, Domain* cd_3, Domain* cd_4) {

        VectorThree ndr_1 {cd_2->m_pos - cd_1->m_pos};
        VectorThree ndr_2 {cd_3->m_pos - cd_2->m_pos};
        VectorThree ndr_3 {cd_4->m_pos - cd_3->m_pos};

        // Check that the domains are not on opposite sides of the junction
        bool junction_ok {true};
        if ((ndr_1 != ndr_2 and ndr_1 != -ndr_2) and (ndr_1 == -ndr_3)) {
            junction_ok = false;
        }

        return junction_ok;
    }

    bool FlexibleBindingPotential::check_edge_pair_junction(Domain* cd_1,
            Domain* cd_2, int i) {

        Domain* cd_j1;
        Domain* cd_j2;
        Domain* cd_j3;
        Domain* cd_j4;
        if (i == -1) {
            cd_j1 = *cd_1 + -2;
            cd_j2 = *cd_1 + -1;
            cd_j3 = cd_1;
            cd_j4 = cd_2;
        }
        else {
            cd_j1 = cd_1;
            cd_j2 = cd_2;
            cd_j3 = *cd_2 + 1;
            cd_j4 = *cd_2 + 2;
        }
        bool junction_ok {true};
        if (check_domains_exist_and_bound({cd_j1, cd_j2, cd_j3, cd_j4})) {
            Domain& cd_bound_j2 {*cd_j2->m_bound_domain};
            Domain& cd_bound_j3 {*cd_j3->m_bound_domain};
            bool bound_same_chain {cd_bound_j2.m_c == cd_bound_j3.m_c};
            if (bound_same_chain) {
                if ((cd_bound_j2.m_d == cd_bound_j3.m_d + 1)) {
                    if (not check_doubly_contig_junction(cd_j1, cd_j2, cd_j3,
                            cd_j4)) {
                        junction_ok = false;
                    }
                }
            }
        }

        return junction_ok;
    }

    pair<double, int> FlexibleBindingPotential::check_pair_stacking(Domain* cd_1,
            Domain* cd_2) {

        double delta_e {0};
        int stacked {0};
        VectorThree ndr {cd_2->m_pos - cd_1->m_pos};
        if (ndr == cd_1->m_ore or ndr == -cd_1->m_ore) {
            ;
        }
        else if (cd_1->check_twist_constraint(ndr, *cd_2)) {
            delta_e += m_pot.stacking_energy(*cd_1, *cd_2);
            stacked = 1;
        }

        return {delta_e, stacked};
    }

    MisbindingPotential::MisbindingPotential(OrigamiPotential& pot):
            m_pot {pot} {
    }

    double RestrictiveMisbindingPotential::bind_domains(
            Domain& cd_i,
            Domain& cd_j) {
        m_constraints_violated = true;
        if (not check_domain_orientations_opposing(cd_i, cd_j)) {
            return 0;
        }
        m_constraints_violated = false;
        return m_pot.hybridization_energy(cd_i, cd_j);
    }

    double DisallowedMisbindingPotential::bind_domains(
            Domain&, Domain&) {
        m_constraints_violated = true;
        return 0;
    }

    OrigamiPotential::OrigamiPotential(
            const vector<vector<int>> identities,
            const vector<vector<string>>& sequences,
            InputParameters& params) :
            m_energy_filebase {params.m_energy_filebase},
            m_temp {params.m_temp},
            m_cation_M {params.m_cation_M},
            m_identities {identities},
            m_sequences {sequences},
            m_stacking_pot {params.m_stacking_pot} {

        if (params.m_binding_pot == "Restrictive") {
            m_binding_pot = new RestrictiveBindingPotential(*this);
        }
        else if (params.m_binding_pot == "Flexible") {
            m_binding_pot = new FlexibleBindingPotential(*this);
        }
        else {
            std::cout << "No such binding potential";
        }

        if (params.m_misbinding_pot == "Restrictive") {
            m_misbinding_pot = new RestrictiveMisbindingPotential(*this);
        }
        else if (params.m_misbinding_pot == "Disallowed") {
            m_misbinding_pot = new DisallowedMisbindingPotential(*this);
        }
        else {
            std::cout << "No such binding potential";
        }

        if (m_stacking_pot == "Constant") {
            m_stacking_ene = params.m_stacking_ene;
        }

        get_energies();
    }

    OrigamiPotential::~OrigamiPotential() {
        delete m_binding_pot;
        delete m_misbinding_pot;
    }

    void OrigamiPotential::update_temp(double temp) {

        // Update hybridization and stacking energy tables
        m_temp = temp;
        if (m_hybridization_energy_tables.count(temp) == 0) {
            get_energies();
            m_hybridization_energy_tables[temp] = m_hybridization_energies;
            m_hybridization_enthalpy_tables[temp] = m_hybridization_enthalpies;
            m_hybridization_entropy_tables[temp] = m_hybridization_entropies;
            m_stacking_energy_tables[temp] = m_stacking_energies;
        }
        else {
            m_hybridization_energies = m_hybridization_energy_tables[temp];
            m_hybridization_enthalpies = m_hybridization_enthalpy_tables[temp];
            m_hybridization_entropies = m_hybridization_entropy_tables[temp];
            m_stacking_energies = m_stacking_energy_tables[temp];
        }
    }

    void OrigamiPotential::get_energies() {
        // Get S, H, and G for all possible interactions and store
        /*if (m_energy_filebase.size() != 0) {
            if (not read_energies_from_file()) {
                calc_energies();
            }
        }
        else {
        */
        calc_energies();
        //}
    }

    void OrigamiPotential::calc_energies() {
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

    void OrigamiPotential::calc_energy(string seq_i, string seq_j,
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
        if (m_stacking_pot == "SequenceSpecific") {
            s_energy = nearestNeighbour::calc_seq_spec_stacking_energy(seq_i,
                    seq_j, m_temp, m_cation_M);
        }
        else if (m_stacking_pot == "Constant") {
            s_energy = m_stacking_ene;
        }
        else {
            std::cout << "No such stacking potential";
        }
        m_stacking_energies[key] = s_energy;
    }

    bool OrigamiPotential::read_energies_from_file() {

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

    void OrigamiPotential::write_energies_to_file() {
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

    pair<double, int> OrigamiPotential::bind_domain(Domain& cd_i) {
        m_constraints_violated = false;
        Domain& cd_j {*cd_i.m_bound_domain};
        double delta_e {0};
        int stacked_pairs {0};
        bool comp {check_domains_complementary(cd_i, cd_j)};
        if (comp) {
            std::tie(delta_e, stacked_pairs) = m_binding_pot->bind_domains(
                    cd_i, cd_j);
            m_constraints_violated = m_binding_pot->m_constraints_violated;
        }
        else {
            delta_e += m_misbinding_pot->bind_domains(cd_i, cd_j);
            m_constraints_violated = m_misbinding_pot->m_constraints_violated;
        }

        return {delta_e, stacked_pairs};
    }

    pair<double, int> OrigamiPotential::check_stacking(Domain& cd_i,
            Domain& cd_j) {

        return m_binding_pot->check_stacking(cd_i, cd_j);
    }

    double OrigamiPotential::hybridization_energy(const Domain& cd_i,
            const Domain& cd_j) const {
        pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
        return m_hybridization_energies.at(key);
    }

    double OrigamiPotential::hybridization_enthalpy(const Domain& cd_i,
            const Domain& cd_j) const {
        pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
        return m_hybridization_enthalpies.at(key);
    }

    double OrigamiPotential::hybridization_entropy(const Domain& cd_i,
            const Domain& cd_j) const {
        pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
        return m_hybridization_entropies.at(key);
    }

    double OrigamiPotential::stacking_energy(const Domain& cd_i, const Domain& cd_j) const {
        pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
        return m_stacking_energies.at(key);
    }

    bool OrigamiPotential::check_domains_complementary(Domain& cd_i, Domain& cd_j) {
        bool comp;
        if (cd_i.m_d_ident == -cd_j.m_d_ident) {
            comp = true;
        }
        else {
            comp = false;
        }
        return comp;
    }
}
