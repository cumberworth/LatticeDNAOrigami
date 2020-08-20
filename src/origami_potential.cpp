// origami_potential.cpp

#include "origami_potential.h"
#include "nearest_neighbour.h"

#include <fstream>
#include <iostream>
#include <utility>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/utility.hpp>

namespace potential {

using std::cout;

using nearestNeighbour::calc_comp_seq;
using utility::Occupancy;
using utility::OrigamiMisuse;

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

bool doubly_contiguous(Domain* cd_1, Domain* cd_2) {

    if (cd_1->m_c != cd_2->m_c or cd_2->m_d != cd_1->m_d + 1) {
        return false;
    }

    Domain* cd_bound_1 {cd_1->m_bound_domain};
    Domain* cd_bound_2 {cd_2->m_bound_domain};
    if (cd_bound_1->m_c != cd_bound_2->m_c) {
        return false;
    }

    if (abs(cd_bound_2->m_d - cd_bound_1->m_d) == 1) {
        return true;
    }

    return false;
}

bool check_pair_stacked(Domain* cd_1, Domain* cd_2) {
    bool stacked {false};
    if (cd_1->m_d > cd_2->m_d) {
        Domain* hold {cd_1};
        cd_1 = cd_2;
        cd_2 = hold;
    }
    VectorThree ndr {cd_2->m_pos - cd_1->m_pos};
    if (ndr != cd_1->m_ore and ndr != -cd_1->m_ore) {
        stacked = cd_1->check_twist_constraint(ndr, *cd_2);
    }

    return stacked;
}

int check_junction_stacking_penalty(
        Domain& cd_j1,
        Domain& cd_j2,
        Domain& cd_j3,
        Domain& cd_j4,
        Domain& cd_k1,
        Domain& cd_k2) {

    int stacking_penalty {0};
    VectorThree ndr_k1 {cd_k2.m_pos - cd_k1.m_pos};
    if (ndr_k1 == cd_k1.m_ore) {
        VectorThree ndr_1 {cd_j2.m_pos - cd_j1.m_pos};
        VectorThree ndr_3 {cd_j4.m_pos - cd_j3.m_pos};
        if (ndr_1 == ndr_3) {
            stacking_penalty = 2;
        }
        else if (ndr_1 == -ndr_3) {
            stacking_penalty = 0;
        }
        else {
            stacking_penalty = 1;
        }
        //        cout << stacking_penalty << " (" << cd_j1.m_c << " " <<
        //        cd_j1.m_d
        //             << "), (" << cd_j2.m_c << " " << cd_j2.m_d << "), (" <<
        //             cd_j3.m_c
        //             << " " << cd_j3.m_d << "), (" << cd_j4.m_c << " " <<
        //             cd_j4.m_d
        //             << "), (" << cd_k1.m_c << " " << cd_k1.m_d << "), (" <<
        //             cd_k2.m_c
        //             << " " << cd_k2.m_d << ")\n";
    }
    return stacking_penalty;
}

BindingPotential::BindingPotential(OrigamiPotential& pot): m_pot {pot} {}

DeltaConfig BindingPotential::bind_domains(Domain& cd_i, Domain& cd_j) {

    m_delta_config = {};
    m_constraints_violated = false;
    if (not check_domain_orientations_opposing(cd_i, cd_j)) {
        m_constraints_violated = true;
        return {};
    }

    calc_stacking_and_steric_terms(cd_i, cd_j);
    if (m_constraints_violated) {
        return {};
    }
    m_delta_config.e += m_pot.hybridization_energy(cd_i, cd_j);

    return m_delta_config;
}

DeltaConfig BindingPotential::check_stacking(Domain& cd_i, Domain& cd_j) {

    m_delta_config = {};
    m_constraints_violated = false;
    calc_stacking_and_steric_terms(cd_i, cd_j);

    return m_delta_config;
}

void BindingPotential::check_triplet_single_stacking(
        Domain* cd_h1,
        Domain* cd_h2,
        Domain* cd_h3) {

    VectorThree ndr_1 {cd_h2->m_pos - cd_h1->m_pos};
    if (ndr_1 == cd_h1->m_ore) {
        return;
    }
    VectorThree ndr_2 {cd_h3->m_pos - cd_h2->m_pos};
    if (ndr_1 != ndr_2) {
        m_delta_config.e -= m_pot.stacking_energy(*cd_h2, *cd_h3);
        m_delta_config.stacked_pairs -= 1;
        //        cout << "(" << cd_h1->m_c << " " << cd_h1->m_d << "), (" <<
        //        cd_h2->m_c
        //             << " " << cd_h2->m_d << "), (" << cd_h3->m_c << " " <<
        //             cd_h3->m_d
        //             << ")\n";
    }
}

void BindingPotential::check_triplet_double_stacking(
        Domain* cd_h1,
        Domain* cd_h2,
        Domain* cd_h3) {

    VectorThree ndr_1 {cd_h2->m_pos - cd_h1->m_pos};
    VectorThree ndr_2 {cd_h3->m_pos - cd_h2->m_pos};
    if (ndr_1 != ndr_2) {
        m_delta_config.e -= m_pot.stacking_energy(*cd_h1, *cd_h2) / 2;
        m_delta_config.e -= m_pot.stacking_energy(*cd_h2, *cd_h3) / 2;
        m_delta_config.stacked_pairs -= 1;
        //        cout << "(" << cd_h1->m_c << " " << cd_h1->m_d << "), (" <<
        //        cd_h2->m_c
        //             << " " << cd_h2->m_d << "), (" << cd_h3->m_c << " " <<
        //             cd_h3->m_d
        //             << ")\n";
    }
}

void BindingPotential::check_triply_contig_helix(
        Domain* cd_h1,
        Domain* cd_h2,
        Domain* cd_h3) {

    VectorThree ndr_1 {cd_h2->m_pos - cd_h1->m_pos};
    VectorThree ndr_2 {cd_h3->m_pos - cd_h2->m_pos};
    if (ndr_1 != ndr_2) {
        m_constraints_violated = true;
    }
}

void JunctionBindingPotential::calc_stacking_and_steric_terms(
        Domain& cd_i,
        Domain& cd_j) {

    int j {0};
    for (auto cd: {&cd_i, &cd_j}) {
        check_constraints(cd, j);
        if (m_constraints_violated) {
            return;
        }
        j++;
    }
    check_central_triplet_stacking_combos(cd_i, cd_j);
}

void JunctionBindingPotential::check_constraints(Domain* cd, int j) {
    for (int i: {-1, 0}) {
        Domain* cd_1 {*cd + i};
        Domain* cd_2 {*cd + (i + 1)};
        if (not check_domains_exist_and_bound({cd_1, cd_2})) {
            continue;
        }
        Domain& cd_bound_1 {*cd_1->m_bound_domain};
        Domain& cd_bound_2 {*cd_2->m_bound_domain};
        bool bound_same_chain {cd_bound_1.m_c == cd_bound_2.m_c};
        if (bound_same_chain) {

            // Same helix case
            if (cd_bound_1.m_d == cd_bound_2.m_d - 1) {
                check_doubly_contig_helix_pair(cd_1, cd_2, i, j);
            }

            // Crossover case
            else if (cd_bound_1.m_d == cd_bound_2.m_d + 1) {
                check_doubly_contig_junction_pair(cd_1, cd_2);
            }

            // Not doubly contiguous
            else {
                check_regular_pair_constraints(cd_1, cd_2, i);
            }
        }
        else {
            check_regular_pair_constraints(cd_1, cd_2, i);
        }
        if (m_constraints_violated) {
            return;
        }
    }
    Domain* cd_prev {*cd + -1};
    Domain* cd_forw {*cd + 1};
    if (check_domains_exist_and_bound({cd_prev, cd_forw})) {
        if (check_pair_stacked(cd, cd_forw)) {
            if (check_pair_stacked(cd_prev, cd)) {
                if (doubly_contiguous(cd_prev, cd) and
                    doubly_contiguous(cd, cd_forw)) {
                    check_triply_contig_helix(cd_prev, cd, cd_forw);
                    if (m_constraints_violated) {
                        return;
                    }
                }
                check_triplet_double_stacking(cd_prev, cd, cd_forw);
            }
            else {
                check_triplet_single_stacking(cd_prev, cd, cd_forw);
            }
        }
    }
}

void JunctionBindingPotential::check_regular_pair_constraints(
        Domain* cd_1,
        Domain* cd_2,
        int i) {

    VectorThree ndr {cd_2->m_pos - cd_1->m_pos};
    if (not cd_1->check_kink_constraint(ndr, *cd_2)) {
        m_constraints_violated = true;
        return;
    }
    if (check_pair_stacked(cd_1, cd_2)) {
        m_delta_config.e += m_pot.stacking_energy(*cd_1, *cd_2);
        m_delta_config.stacked_pairs += 1;
        //        cout << cd_1->m_c << " " << cd_1->m_d << ", " << cd_2->m_c <<
        //        " "
        //             << cd_2->m_d << "\n";
        if (i == -1) {
            check_backward_single_junction(cd_1, cd_2);
        }
        else {
            check_forward_single_junction(cd_1, cd_2);
        }
    }
    else {
        check_central_single_junction(cd_1, cd_2);
    }

    if (i == -1) {
        check_backward_triplet_stacking_combos(cd_1, cd_2, i);
    }
    else {
        check_forward_triplet_stacking_combos(cd_1, cd_2, i);
    }
}

void JunctionBindingPotential::check_doubly_contig_helix_pair(
        Domain* cd_1,
        Domain* cd_2,
        int i,
        int j) {

    VectorThree ndr {cd_2->m_pos - cd_1->m_pos};
    if (ndr == cd_1->m_ore or ndr == -cd_1->m_ore) {
        m_constraints_violated = true;
    }
    else if (cd_1->check_twist_constraint(ndr, *cd_2)) {
        if (j == 0) {
            m_delta_config.e += m_pot.stacking_energy(*cd_1, *cd_2);
            m_delta_config.stacked_pairs += 1;
            //            cout << cd_1->m_c << " " << cd_1->m_d << ", " <<
            //            cd_2->m_c << " "
            //                 << cd_2->m_d << "\n";
        }
        if (i == -1) {
            check_backward_triplet_stacking_combos(cd_1, cd_2, i);
            check_backward_single_junction(cd_1, cd_2);
        }
        else {
            check_forward_triplet_stacking_combos(cd_1, cd_2, i);
            check_forward_single_junction(cd_1, cd_2);
        }
    }
    else {
        m_constraints_violated = true;
    }
}

void JunctionBindingPotential::check_doubly_contig_junction_pair(
        Domain* cd_1,
        Domain* cd_2) {

    Domain* cd_j1 {*cd_1 + -1};
    Domain* cd_j2 {cd_1};
    Domain* cd_j3;
    Domain* cd_j4;
    Domain* cd_k1 {cd_1};
    Domain* cd_k2 {cd_2};
    VectorThree ndr {cd_k2->m_pos - cd_k1->m_pos};
    if (cd_k1->m_ore != ndr) {
        m_constraints_violated = true;
        return;
    }
    cd_j3 = cd_k2;
    cd_j4 = *cd_k2 + 1;
    if (check_domains_exist_and_bound({cd_j4})) {
        if (check_domains_exist_and_bound({cd_j1})) {
            check_junction(cd_j1, cd_j2, cd_j3, cd_j4, cd_k1, cd_k2);
        }
        if (check_pair_stacked(cd_j3, cd_j4)) {
            check_triplet_single_stacking(cd_k1, cd_k2, cd_j4);
        }
    }

    cd_j3 = cd_k2->m_bound_domain;
    cd_j4 = *cd_j3 + -1;
    if (check_domains_exist_and_bound({cd_j4})) {
        if (check_domains_exist_and_bound({cd_j1})) {
            check_junction(cd_j1, cd_j2, cd_j3, cd_j4, cd_k1, cd_k2);
        }
        if (check_pair_stacked(cd_j4, cd_j3)) {
            check_triplet_single_stacking(cd_k1, cd_k2, cd_j4);
        }
    }
}

void JunctionBindingPotential::check_backward_triplet_stacking_combos(
        Domain* cd_1,
        Domain* cd_2,
        int i) {

    Domain* cd_h1;
    Domain* cd_h2;
    Domain* cd_h3;
    cd_h2 = cd_1;
    cd_h3 = cd_2;
    if (not check_pair_stacked(cd_h2, cd_h3)) {
        return;
    }
    bool h2_h3_doubly_contig {doubly_contiguous(cd_h2, cd_h3)};

    // Check same chain
    cd_h1 = *cd_1 + -1;
    if (check_domains_exist_and_bound({cd_h1})) {
        if (check_pair_stacked(cd_h1, cd_h2)) {
            if (h2_h3_doubly_contig and doubly_contiguous(cd_h1, cd_h2)) {
                check_triply_contig_helix(cd_h1, cd_h2, cd_h3);
                if (m_constraints_violated) {
                    return;
                }
            }
            check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
        }
        else {
            check_triplet_single_stacking(cd_h1, cd_h2, cd_h3);
        }
    }

    // Check helices that extend to bound chain
    Domain* cd_h2_prev {cd_h1};
    cd_h1 = *(cd_1->m_bound_domain) + 1;
    if (check_domains_exist_and_bound({cd_h1}) and
        cd_h1->m_bound_domain != cd_h3 and not h2_h3_doubly_contig) {
        if (check_pair_stacked(cd_1->m_bound_domain, cd_h1)) {
            check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
        }
    }
    cd_h1 = *(cd_1->m_bound_domain) + -1;
    if (check_domains_exist_and_bound({cd_h1}) and
        cd_h1->m_bound_domain != cd_h2_prev and not h2_h3_doubly_contig) {
        if (check_pair_stacked(cd_h1, cd_1->m_bound_domain)) {
            check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
        }
        else {
            check_triplet_single_stacking(cd_h1, cd_1->m_bound_domain, cd_h3);
        }
    }
}

void JunctionBindingPotential::check_forward_triplet_stacking_combos(
        Domain* cd_1,
        Domain* cd_2,
        int i) {

    Domain* cd_h1;
    Domain* cd_h2;
    Domain* cd_h3;
    cd_h2 = cd_2;
    cd_h1 = cd_1;
    bool first_pair_stacked {check_pair_stacked(cd_h1, cd_h2)};
    bool h1_h2_doubly_contig {doubly_contiguous(cd_h1, cd_h2)};

    // Check same chain
    cd_h3 = *(cd_2) + 1;
    if (check_domains_exist_and_bound({cd_h3})) {
        bool second_pair_stacked {check_pair_stacked(cd_h2, cd_h3)};
        if (first_pair_stacked and second_pair_stacked) {
            if (h1_h2_doubly_contig and doubly_contiguous(cd_h2, cd_h3)) {
                check_triply_contig_helix(cd_h1, cd_h2, cd_h3);
                if (m_constraints_violated) {
                    return;
                }
            }
            check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
        }
        else if (second_pair_stacked) {
            check_triplet_single_stacking(cd_h1, cd_h2, cd_h3);
        }
    }

    // Check helices that extend to bound chain
    Domain* cd_h2_next {cd_h3};
    cd_h3 = *(cd_2->m_bound_domain) + 1;
    if (check_domains_exist_and_bound({cd_h3}) and
        cd_h3->m_bound_domain != cd_h2_next and not h1_h2_doubly_contig) {
        if (check_pair_stacked(cd_2->m_bound_domain, cd_h3)) {
            if (first_pair_stacked) {
                check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
            }
            else {
                check_triplet_single_stacking(cd_h1, cd_h2, cd_h3);
            }
        }
    }
    cd_h3 = *(cd_2->m_bound_domain) + -1;
    if (check_domains_exist_and_bound({cd_h3}) and
        cd_h3->m_bound_domain != cd_h1 and not h1_h2_doubly_contig) {
        if (check_pair_stacked(cd_h3, cd_2->m_bound_domain)) {
            if (first_pair_stacked) {
                check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
            }
            else {
                check_triplet_single_stacking(cd_h1, cd_h2, cd_h3);
            }
        }
        else if (first_pair_stacked) {
            check_triplet_single_stacking(cd_h3, cd_2->m_bound_domain, cd_h1);
        }
    }
}

void JunctionBindingPotential::check_central_triplet_stacking_combos(
        Domain& cd_i,
        Domain& cd_j) {

    Domain* cd_h1 {cd_i + -1};
    Domain* cd_h2 {&cd_i};
    Domain* cd_h3;
    cd_h3 = cd_j + 1;
    Domain* cd_h2_next {cd_i + 1};
    if (check_domains_exist_and_bound({cd_h1, cd_h3}) and
        cd_h3->m_bound_domain != cd_h2_next and
        (not doubly_contiguous(cd_h1, cd_h2))) {
        if (check_pair_stacked(&cd_j, cd_h3)) {
            if (check_pair_stacked(cd_h1, cd_h2)) {
                check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
            }
            else {
                check_triplet_single_stacking(cd_h1, cd_h2, cd_h3);
            }
        }
    }
    cd_h3 = cd_j + -1;
    if (check_domains_exist_and_bound({cd_h1, cd_h3}) and
        cd_h3->m_bound_domain != cd_h1 and
        (not doubly_contiguous(cd_h1, cd_h2))) {
        if (check_pair_stacked(cd_h1, cd_h2)) {
            if (check_pair_stacked(cd_h3, &cd_j)) {
                check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
            }
            else {
                check_triplet_single_stacking(cd_h3, &cd_j, cd_h1);
            }
        }
        else if (check_pair_stacked(cd_h3, &cd_j)) {
            check_triplet_single_stacking(cd_h1, cd_h2, cd_h3);
        }
    }

    cd_h2 = &cd_i;
    cd_h3 = cd_i + 1;
    cd_h1 = cd_j + 1;
    Domain* cd_h2_prev {cd_i + -1};
    if (check_domains_exist_and_bound({cd_h1, cd_h3}) and
        cd_h1->m_bound_domain != cd_h3 and
        (not doubly_contiguous(cd_h2, cd_h3))) {
        if (check_pair_stacked(&cd_j, cd_h1) and
            check_pair_stacked(cd_h2, cd_h3)) {
            check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
        }
    }
    cd_h1 = cd_j + -1;
    if (check_domains_exist_and_bound({cd_h1, cd_h3}) and
        cd_h1->m_bound_domain != cd_h2_prev and
        (not doubly_contiguous(cd_h2, cd_h3))) {
        if (check_pair_stacked(cd_h2, cd_h3)) {
            if (check_pair_stacked(cd_h1, &cd_j)) {
                check_triplet_double_stacking(cd_h1, cd_h2, cd_h3);
            }
            else {
                check_triplet_single_stacking(cd_h1, &cd_j, cd_h3);
            }
        }
    }
}

void JunctionBindingPotential::check_backward_single_junction(
        Domain* cd_1,
        Domain* cd_2) {

    // Passed domains are the second junction pair
    Domain* cd_j1;
    Domain* cd_j2;
    Domain* cd_j3 {cd_1};
    Domain* cd_j4 {cd_2};
    Domain* cd_k1;
    Domain* cd_k2;

    // Find possible kink pair combinations
    vector<pair<Domain*, Domain*>> first_sel {};

    // Add kink pair that is the same chain as the first junction pair
    Domain* cd_j3_bac {cd_j3->m_backward_domain};
    first_sel.push_back({cd_j3, cd_j3_bac});

    // Add kink pairs that are on the chain bound to j3
    Domain* cd_j3_bound {cd_j3->m_bound_domain};
    Domain* cd_j4_bound {cd_j4->m_bound_domain};
    Domain* cd_j3_bound_for {cd_j3_bound->m_forward_domain};
    Domain* cd_j3_bound_bac {cd_j3_bound->m_backward_domain};
    if (cd_j3->m_bound_domain->m_c == cd_j4->m_bound_domain->m_c) {
        int diff {cd_j4_bound->m_d - cd_j3_bound->m_d};

        // Must prevent double counting
        if (abs(diff) == 1 and cd_j3->m_c > cd_j3_bound->m_c) {
            return;
        }
        if (diff == 1) {
            first_sel.push_back({cd_j3_bound, cd_j3_bound_bac});
        }
        else {
            first_sel.push_back({cd_j3_bound, cd_j3_bound_for});
            first_sel.push_back({cd_j3_bound, cd_j3_bound_bac});
        }
    }
    else if (
            check_domains_exist_and_bound(
                    {cd_j3_bac, cd_j3_bound_bac}) and
            cd_j3_bac->m_bound_domain == cd_j3_bound_bac) {
        first_sel.push_back({cd_j3_bound, cd_j3_bound_for});
    }
    else if (
            check_domains_exist_and_bound(
                    {cd_j3_bac, cd_j3_bound_for}) and
            cd_j3_bac->m_bound_domain == cd_j3_bound_for) {
        first_sel.push_back({cd_j3_bound, cd_j3_bound_bac});
    }
    else {
        first_sel.push_back({cd_j3_bound, cd_j3_bound_for});
        first_sel.push_back({cd_j3_bound, cd_j3_bound_bac});
    }
    for (auto sel1: first_sel) {
        cd_k2 = sel1.first;
        cd_k1 = sel1.second;
        if (not check_domains_exist_and_bound({cd_k1})) {
            continue;
        }

        // If kink is stacked, not a kink
        if (check_pair_stacked(cd_k1, cd_k2)) {
            continue;
        }

        // Find possible second junction pairs
        vector<pair<Domain*, Domain*>> second_sel {};

        // Add junction that is on the same chain as the kink pair
        int dir {cd_k1->m_d - cd_k2->m_d};
        Domain* cd_k1_next {*cd_k1 + dir};
        second_sel.push_back({cd_k1, cd_k1_next});

        // Add junction pairs that are on chain bound to k2
        Domain* cd_k1_bound {cd_k1->m_bound_domain};
        Domain* cd_k1_bound_for {cd_k1_bound->m_forward_domain};
        Domain* cd_k1_bound_bac {cd_k1_bound->m_backward_domain};
        if (check_domains_exist_and_bound({cd_k1_next})) {

            // Prevent double counting
            Domain* cd_k1_next_bound {cd_k1_next->m_bound_domain};
            if (cd_k1_bound_for != cd_k1_next_bound) {
                second_sel.push_back({cd_k1_bound, cd_k1_bound_for});
            }
            if (cd_k1_bound_bac != cd_k1_next_bound) {
                second_sel.push_back({cd_k1_bound, cd_k1_bound_bac});
            }
        }
        else {
            second_sel.push_back({cd_k1_bound, cd_k1_bound_for});
            second_sel.push_back({cd_k1_bound, cd_k1_bound_bac});
        }
        for (auto sel2: second_sel) {
            cd_j2 = sel2.first;
            cd_j1 = sel2.second;
            if (not check_domains_exist_and_bound({cd_j1})) {
                continue;
            }
            check_junction(cd_j1, cd_j2, cd_j3, cd_j4, cd_k1, cd_k2);
        }
    }
}

void JunctionBindingPotential::check_forward_single_junction(
        Domain* cd_1,
        Domain* cd_2) {

    // Passed domains are the first junction pair
    Domain* cd_j1 {cd_1};
    Domain* cd_j2 {cd_2};
    Domain* cd_j3;
    Domain* cd_j4;
    Domain* cd_k1;
    Domain* cd_k2;

    // Find possible kink pair combinations
    vector<pair<Domain*, Domain*>> first_sel {};

    // Add kink pair that is the same chain as the first junction pair
    Domain* cd_j2_for {cd_j2->m_forward_domain};
    first_sel.push_back({cd_j2, cd_j2_for});

    // Add kink pairs that are on the chain bound to j2
    Domain* cd_j1_bound {cd_j1->m_bound_domain};
    Domain* cd_j2_bound {cd_j2->m_bound_domain};
    Domain* cd_j2_bound_for {cd_j2_bound->m_forward_domain};
    Domain* cd_j2_bound_bac {cd_j2_bound->m_backward_domain};
    if (cd_j1_bound->m_c == cd_j2_bound->m_c) {
        int diff {cd_j2_bound->m_d - cd_j1_bound->m_d};

        // Must prevent double counting
        if (abs(diff) == 1 and cd_j1->m_c > cd_j1_bound->m_c) {
            return;
        }
        if (diff == 1) {
            first_sel.push_back({cd_j2_bound, cd_j2_bound_for});
        }
        else {
            first_sel.push_back({cd_j2_bound, cd_j2_bound_for});
            first_sel.push_back({cd_j2_bound, cd_j2_bound_bac});
        }
    }
    // What if it is doubly contig on both sides?
    else if (
            check_domains_exist_and_bound(
                    {cd_j2_for, cd_j2_bound_bac}) and
            cd_j2_for->m_bound_domain == cd_j2_bound_bac) {
        first_sel.push_back({cd_j2_bound, cd_j2_bound_for});
    }
    else if (
            check_domains_exist_and_bound(
                    {cd_j2_for, cd_j2_bound_for}) and
            cd_j2_for->m_bound_domain == cd_j2_bound_for) {
        first_sel.push_back({cd_j2_bound, cd_j2_bound_bac});
    }
    else {
        first_sel.push_back({cd_j2_bound, cd_j2_bound_for});
        first_sel.push_back({cd_j2_bound, cd_j2_bound_bac});
    }
    for (auto sel1: first_sel) {
        cd_k1 = sel1.first;
        cd_k2 = sel1.second;
        if (not check_domains_exist_and_bound({cd_k2})) {
            continue;
        }

        // If kink is stacked, not a kink
        if (check_pair_stacked(cd_k1, cd_k2)) {
            continue;
        }

        // Find possible second junction pairs
        vector<pair<Domain*, Domain*>> second_sel {};

        // Add junction that is on the same chain as the kink pair
        int dir {cd_k2->m_d - cd_k1->m_d};
        Domain* cd_k2_next {*cd_k2 + dir};
        second_sel.push_back({cd_k2, cd_k2_next});

        // Add junction pairs that are on chain bound to k2
        Domain* cd_k2_bound {cd_k2->m_bound_domain};
        Domain* cd_k2_bound_for {cd_k2_bound->m_forward_domain};
        Domain* cd_k2_bound_bac {cd_k2_bound->m_backward_domain};
        if (check_domains_exist_and_bound({cd_k2_next})) {

            // Prevent double counting
            Domain* cd_k2_next_bound {cd_k2_next->m_bound_domain};
            if (cd_k2_bound_for != cd_k2_next_bound) {
                second_sel.push_back({cd_k2_bound, cd_k2_bound_for});
            }
            if (cd_k2_bound_bac != cd_k2_next_bound) {
                second_sel.push_back({cd_k2_bound, cd_k2_bound_bac});
            }
        }
        else {
            second_sel.push_back({cd_k2_bound, cd_k2_bound_for});
            second_sel.push_back({cd_k2_bound, cd_k2_bound_bac});
        }
        for (auto sel2: second_sel) {
            cd_j3 = sel2.first;
            cd_j4 = sel2.second;
            if (not check_domains_exist_and_bound({cd_j4})) {
                continue;
            }
            check_junction(cd_j1, cd_j2, cd_j3, cd_j4, cd_k1, cd_k2);
        }
    }
}

void JunctionBindingPotential::check_central_single_junction(
        Domain* cd_1,
        Domain* cd_2) {

    // Passed domains are the kink pair
    Domain* cd_j1;
    Domain* cd_j2;
    Domain* cd_j3;
    Domain* cd_j4;
    Domain* cd_k1 {cd_1};
    Domain* cd_k2 {cd_2};

    // Kinked pair cannot be doubly contiguous
    if (cd_k1->m_bound_domain->m_c == cd_k2->m_bound_domain->m_c and
        abs(cd_k1->m_bound_domain->m_d - cd_k2->m_bound_domain->m_d) == 1) {
        return;
    }

    // Find possible first junction pairs
    vector<pair<Domain*, Domain*>> first_sel {};

    // Add first junction pair that is on the same chain as the kink pair
    Domain* cd_k1_bac {cd_k1->m_backward_domain};
    first_sel.push_back({cd_k1, cd_k1_bac});

    // Add first junction pairs bound to k1
    Domain* cd_k1_bound {cd_k1->m_bound_domain};
    Domain* cd_k1_bound_for {cd_k1_bound->m_forward_domain};
    Domain* cd_k1_bound_bac {cd_k1_bound->m_backward_domain};
    if (check_domains_exist_and_bound({cd_k1_bac})) {

        // If j1 and j2 are doubly contiguous, then this will be double
        // counted
        Domain* cd_k1_bac_bound {cd_k1_bac->m_bound_domain};
        if (cd_k1_bound_for != cd_k1_bac_bound) {
            first_sel.push_back({cd_k1_bound, cd_k1_bound_for});
        }
        if (cd_k1_bound_bac != cd_k1_bac_bound) {
            first_sel.push_back({cd_k1_bound, cd_k1_bound_bac});
        }
    }
    else {
        first_sel.push_back({cd_k1_bound, cd_k1_bound_for});
        first_sel.push_back({cd_k1_bound, cd_k1_bound_bac});
    }
    for (auto sel1: first_sel) {
        cd_j2 = sel1.first;
        cd_j1 = sel1.second;
        if (not check_domains_exist_and_bound({cd_j1})) {
            continue;
        }

        // Find possible second junction pairs
        vector<pair<Domain*, Domain*>> second_sel {};

        // Add junction pair that is on the same chain as the kink pair
        Domain* cd_k2_for {cd_k2->m_forward_domain};
        second_sel.push_back({cd_k2, cd_k2_for});

        // Add junction pairs that are on chain bound to k2
        Domain* cd_k2_bound {cd_k2->m_bound_domain};
        Domain* cd_k2_bound_for {cd_k2_bound->m_forward_domain};
        Domain* cd_k2_bound_bac {cd_k2_bound->m_backward_domain};
        if (check_domains_exist_and_bound({cd_k2_for})) {

            // If j3 and j4 are doubly contiguous, then this will be double
            // counted
            Domain* cd_k2_for_bound {cd_k2_for->m_bound_domain};
            if (cd_k2_bound_for != cd_k2_for_bound) {
                second_sel.push_back({cd_k2_bound, cd_k2_bound_for});
            }
            if (cd_k2_bound_bac != cd_k2_for_bound) {
                second_sel.push_back({cd_k2_bound, cd_k2_bound_bac});
            }
        }
        else {
            second_sel.push_back({cd_k2_bound, cd_k2_bound_for});
            second_sel.push_back({cd_k2_bound, cd_k2_bound_bac});
        }
        for (auto sel2: second_sel) {
            cd_j3 = sel2.first;
            cd_j4 = sel2.second;
            if (not check_domains_exist_and_bound({cd_j4})) {
                continue;
            }
            check_junction(cd_j1, cd_j2, cd_j3, cd_j4, cd_k1, cd_k2);
        }
    }
}

void JunctionBindingPotential::check_junction(
        Domain* cd_j1,
        Domain* cd_j2,
        Domain* cd_j3,
        Domain* cd_j4,
        Domain* cd_k1,
        Domain* cd_k2) {

    // Put domains in order along chain
    Domain* hold;
    if (cd_k1->m_d > cd_k2->m_d) {
        hold = cd_k1;
        cd_k1 = cd_k2;
        cd_k2 = hold;
        hold = cd_j1;
        cd_j1 = cd_j4;
        cd_j4 = hold;
        hold = cd_j2;
        cd_j2 = cd_j3;
        cd_j3 = hold;
    }
    // this is inefficient for three quarter domains
    if (not(check_pair_stacked(cd_j1, cd_j2) and
            check_pair_stacked(cd_j3, cd_j4))) {
        return;
    }
    if (not cd_j1->check_junction_constraint(
                *cd_j2, *cd_j3, *cd_j4, *cd_k1, *cd_k2)) {
        m_constraints_violated = true;
        return;
    }
    auto stacking_penalty {check_junction_stacking_penalty(
            *cd_j1, *cd_j2, *cd_j3, *cd_j4, *cd_k1, *cd_k2)};
    if (stacking_penalty == 1) {
        m_delta_config.e -= m_pot.stacking_energy(*cd_j1, *cd_j2) / 2;
        m_delta_config.e -= m_pot.stacking_energy(*cd_j3, *cd_j4) / 2;
        m_delta_config.stacked_pairs -= 1;
    }
    else if (stacking_penalty == 2) {
        m_delta_config.e -= m_pot.stacking_energy(*cd_j1, *cd_j2);
        m_delta_config.e -= m_pot.stacking_energy(*cd_j3, *cd_j4);
        m_delta_config.stacked_pairs -= 2;
    }
}

MisbindingPotential::MisbindingPotential(OrigamiPotential& pot): m_pot {pot} {}

double OpposingMisbindingPotential::bind_domains(Domain& cd_i, Domain& cd_j) {
    m_constraints_violated = true;
    if (not check_domain_orientations_opposing(cd_i, cd_j)) {
        return 0;
    }
    m_constraints_violated = false;
    return m_pot.hybridization_energy(cd_i, cd_j);
}

double DisallowedMisbindingPotential::bind_domains(Domain&, Domain&) {
    m_constraints_violated = true;
    return 0;
}

OrigamiPotential::OrigamiPotential(
        const vector<vector<int>> identities,
        const vector<vector<string>>& sequences,
        const vector<double> enthalpies,
        const vector<double> entropies,
        InputParameters& params):
        m_energy_filebase {params.m_energy_filebase},
        m_temp {params.m_temp},
        m_cation_M {params.m_cation_M},
        m_identities {identities},
        m_sequences {sequences},
        m_complementary_enthalpies {enthalpies},
        m_complementary_entropies {entropies},
        m_stacking_pot {params.m_stacking_pot},
        m_hybridization_pot {params.m_hybridization_pot},
        m_apply_mean_field_cor {params.m_apply_mean_field_cor} {

    if (params.m_binding_pot == "FourBody") {
        m_binding_pot = new JunctionBindingPotential(*this);
    }
    else {
        std::cout << "No such binding potential";
    }

    if (params.m_misbinding_pot == "Opposing") {
        m_misbinding_pot = new OpposingMisbindingPotential(*this);
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

    if (m_hybridization_pot == "Uniform") {
        m_binding_h = params.m_binding_h;
        m_binding_s = params.m_binding_s;
        m_misbinding_h = params.m_misbinding_h;
        m_misbinding_s = params.m_misbinding_s;
    }
    if (m_hybridization_pot == "Specified") {
        m_misbinding_h = params.m_misbinding_h;
        m_misbinding_s = params.m_misbinding_s;
    }

    get_energies();
}

OrigamiPotential::~OrigamiPotential() {
    delete m_binding_pot;
    delete m_misbinding_pot;
}

void OrigamiPotential::update_temp(double temp, double stacking_mult) {

    // Update hybridization and stacking energy tables
    m_temp = temp;
    pair<double, double> key {temp, stacking_mult};
    if (m_stacking_energy_tables.count(key) == 0) {

        // THIS ONLY WORKS FOR CONSTANT STACKING
        double old_stacking_ene {m_stacking_ene};
        m_stacking_ene *= stacking_mult;
        get_energies();
        m_stacking_ene = old_stacking_ene;
        m_hybridization_energy_tables[temp] = m_hybridization_energies;
        m_hybridization_enthalpy_tables[temp] = m_hybridization_enthalpies;
        m_hybridization_entropy_tables[temp] = m_hybridization_entropies;
        m_stacking_energy_tables[key] = m_stacking_energies;
    }
    else {
        m_hybridization_energies = m_hybridization_energy_tables[temp];
        m_hybridization_enthalpies = m_hybridization_enthalpy_tables[temp];
        m_hybridization_entropies = m_hybridization_entropy_tables[temp];
        m_stacking_energies = m_stacking_energy_tables[key];
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

    ThermoOfHybrid DH_DS {nearestNeighbour::calc_unitless_init_thermo(m_temp)};
    m_init_enthalpy = DH_DS.enthalpy;
    m_init_entropy = DH_DS.entropy;
    m_init_energy = m_init_enthalpy - m_init_entropy;

    // Loop through all pairs of sequences
    for (size_t c_i {0}; c_i != m_identities.size(); c_i++) {
        for (size_t c_j {0}; c_j != m_identities.size(); c_j++) {
            size_t c_i_length {m_identities[c_i].size()};
            size_t c_j_length {m_identities[c_j].size()};
            for (size_t d_i {0}; d_i != c_i_length; d_i++) {
                int d_i_ident {m_identities[c_i][d_i]};
                for (size_t d_j {0}; d_j != c_j_length; d_j++) {
                    int d_j_ident {m_identities[c_j][d_j]};
                    pair<int, int> key {d_i_ident, d_j_ident};
                    if (m_hybridization_pot == "NearestNeighbour") {
                        string seq_i {m_sequences[c_i][d_i]};
                        string seq_j {m_sequences[c_j][d_j]};
                        calc_hybridization_energy(seq_i, seq_j, key);
                    }
                    else if (m_hybridization_pot == "Uniform") {
                        calc_hybridization_energy(key);
                    }
                    else if (m_hybridization_pot == "Specified") {
                        set_hybridization_energy(key);
                    }
                    else {
                        std::cout << "No such hybridization potential";
                    }

                    if (m_stacking_pot == "SequenceSpecific") {
                        string seq_i {m_sequences[c_i][d_i]};
                        string seq_j {m_sequences[c_j][d_j]};
                        calc_stacking_energy(seq_i, seq_j, key);
                    }
                    else if (m_stacking_pot == "Constant") {
                        calc_stacking_energy(key);
                    }
                    else {
                        std::cout << "No such stacking potential";
                    }
                }
            }
        }
    }
    write_energies_to_file();
}

void OrigamiPotential::calc_hybridization_energy(
        string seq_i,
        string seq_j,
        pair<int, int> key) {

    double H_hyb {0};
    double S_hyb {0};
    vector<string> comp_seqs {
            nearestNeighbour::find_longest_contig_complement(seq_i, seq_j)};
    int N {0};

    // No interaction if no complementary sequence
    if (comp_seqs.size() == 0) {
        H_hyb = 0;
        S_hyb = 0;
    }
    else {

        // Take average value of H and S of all equal length comp seqs
        for (auto comp_seq: comp_seqs) {
            ThermoOfHybrid DH_DS {
                    nearestNeighbour::calc_unitless_hybridization_thermo(
                            comp_seq, m_temp, m_cation_M)};
            H_hyb += DH_DS.enthalpy - m_init_enthalpy;
            S_hyb += DH_DS.entropy - m_init_entropy;
            N++;
        }
        H_hyb /= N;
        S_hyb /= N;

        S_hyb += log(6);

        if (key.first == -key.second) {
            if (m_apply_mean_field_cor) {
                S_hyb += 3 * log(6);
            }

            // Check that sequences are complementary if they should be
            if (comp_seqs[0].size() != seq_i.size() and
                seq_i.size() == seq_j.size()) {
                cout << "Sequences that should be complementary are not\n";
                throw OrigamiMisuse {};
            }
        }

        // Check that sequences are not complementary if they shouldn't be
        else {
            if (comp_seqs[0].size() == seq_i.size() and
                seq_i.size() == seq_j.size()) {
                cout << "Sequences that should not be complementary are\n";
                throw OrigamiMisuse {};
            }
        }
    }

    m_hybridization_enthalpies[key] = H_hyb;
    m_hybridization_entropies[key] = S_hyb;
    m_hybridization_energies[key] = H_hyb - S_hyb;
}

void OrigamiPotential::calc_hybridization_energy(pair<int, int> key) {
    double H_hyb {0};
    double S_hyb {0};
    if (key.first == -key.second) {
        H_hyb = m_binding_h / m_temp;
        S_hyb = m_binding_s;
        if (m_apply_mean_field_cor) {
            S_hyb += 3 * log(6);
        }
    }
    else {
        H_hyb = m_misbinding_h / m_temp;
        S_hyb = m_misbinding_s;
    }

    S_hyb += log(6);

    m_hybridization_enthalpies[key] = H_hyb;
    m_hybridization_entropies[key] = S_hyb;
    m_hybridization_energies[key] = H_hyb - S_hyb;
}

void OrigamiPotential::set_hybridization_energy(pair<int, int> key) {
    double H_hyb {0};
    double S_hyb {0};
    if (key.first == -key.second) {
        int i {std::abs(key.first) - 1};
        H_hyb = m_complementary_enthalpies[i] / m_temp;
        S_hyb = m_complementary_entropies[i];
        if (m_apply_mean_field_cor) {
            S_hyb += 3 * log(6);
        }
    }
    else {
        H_hyb = m_misbinding_h / m_temp;
        S_hyb = m_misbinding_s;
    }

    S_hyb += log(6);

    m_hybridization_enthalpies[key] = H_hyb;
    m_hybridization_entropies[key] = S_hyb;
    m_hybridization_energies[key] = H_hyb - S_hyb;
}

void OrigamiPotential::calc_stacking_energy(
        string seq_i,
        string seq_j,
        pair<int, int> key) {

    double s_energy {nearestNeighbour::calc_seq_spec_stacking_energy(
            seq_i, seq_j, m_temp, m_cation_M)};
    m_stacking_energies[key] = s_energy;
}

void OrigamiPotential::calc_stacking_energy(pair<int, int> key) {
    m_stacking_energies[key] = m_stacking_ene / m_temp;
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

    bool files_present {
            henergy_file and hhenergy_file and hsenergy_file and senergy_file};
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

DeltaConfig OrigamiPotential::bind_domain(Domain& cd_i) {
    m_constraints_violated = false;
    Domain& cd_j {*cd_i.m_bound_domain};
    bool comp {check_domains_complementary(cd_i, cd_j)};
    DeltaConfig delta_config {};
    if (comp) {
        delta_config = m_binding_pot->bind_domains(cd_i, cd_j);
        m_constraints_violated = m_binding_pot->m_constraints_violated;
    }
    else {
        delta_config.e += m_misbinding_pot->bind_domains(cd_i, cd_j);
        m_constraints_violated = m_misbinding_pot->m_constraints_violated;
    }
    // DEBUG
    /*if (not m_constraints_violated) {
        if (cd_i.m_backward_domain != nullptr and
    cd_i.m_backward_domain->m_state != Occupancy::unassigned) { VectorThree ndr
    {cd_i.m_pos - cd_i.m_backward_domain->m_pos}; if (ndr.abssum() != 1) { cout
    << "Contiguous domains not on adjacent sites\n"; throw OrigamiMisuse {};
            }
        }
        if (cd_i.m_forward_domain != nullptr and cd_i.m_forward_domain->m_state
    != Occupancy::unassigned) { VectorThree ndr {cd_i.m_pos -
    cd_i.m_forward_domain->m_pos}; if (ndr.abssum() != 1) { cout << "Contiguous
    domains not on adjacent sites\n"; throw OrigamiMisuse {};
            }
        }
    }*/

    return delta_config;
}

DeltaConfig OrigamiPotential::check_stacking(Domain& cd_i, Domain& cd_j) {

    return m_binding_pot->check_stacking(cd_i, cd_j);
}

double OrigamiPotential::hybridization_energy(
        const Domain& cd_i,
        const Domain& cd_j) const {
    pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
    return m_hybridization_energies.at(key);
}

double OrigamiPotential::hybridization_enthalpy(
        const Domain& cd_i,
        const Domain& cd_j) const {
    pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
    return m_hybridization_enthalpies.at(key);
}

double OrigamiPotential::hybridization_entropy(
        const Domain& cd_i,
        const Domain& cd_j) const {
    pair<int, int> key {cd_i.m_d_ident, cd_j.m_d_ident};
    return m_hybridization_entropies.at(key);
}

double OrigamiPotential::init_enthalpy() const { return m_init_enthalpy; }

double OrigamiPotential::init_entropy() const { return m_init_entropy; }

double OrigamiPotential::init_energy() const { return m_init_energy; }

double OrigamiPotential::stacking_energy(const Domain& cd_i, const Domain& cd_j)
        const {
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
} // namespace potential
