// nearest_neighbour.h

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "nearest_neighbour.h"

namespace nearestNeighbour {

using std::cout;
using std::log;
using std::reverse;
using std::string;
using std::tie;
using std::tuple;
using std::vector;

double calc_seq_spec_stacking_energy(
        string seq_i,
        string seq_j,
        double,
        double) {

    string nuc_i_back {seq_i.back()};
    string nuc_j_front {seq_j.front()};
    string stack_pair {nuc_i_back + nuc_j_front};
    string comp_stack_pair {calc_comp_seq(stack_pair)};
    string key {stack_pair + "/" + comp_stack_pair};

    return Stacking_Energy.at(key);
}

ThermoOfHybrid calc_unitless_hybridization_thermo(
        string seq,
        double temp,
        double cation_M) {
    ThermoOfHybrid DH_DS {calc_hybridization_H_and_S(seq, cation_M)};

    // Convert from kcal/mol to unitless dimension
    DH_DS.enthalpy = DH_DS.enthalpy * J_Per_Cal * 1000 / R / temp;
    DH_DS.entropy = DH_DS.entropy * J_Per_Cal * 1000 / R;

    return DH_DS;
}

double calc_unitless_hybridization_energy(
        string seq,
        double temp,
        double cation_M) {
    ThermoOfHybrid DH_DS {
            calc_unitless_hybridization_thermo(seq, temp, cation_M)};

    return DH_DS.enthalpy - DH_DS.entropy;
}

ThermoOfHybrid calc_hybridization_H_and_S(string seq, double cation_M) {
    string comp_seq {calc_comp_seq(seq)};

    // Initiation free energy
    double DH_init {NN_Enthalpy.at("INITIATION")};
    double DS_init {NN_Entropy.at("INITIATION")};

    // Symmetry penalty for palindromic sequences
    double DS_sym {0};
    if (seq_is_palindromic(seq)) {
        DS_sym = NN_Entropy.at("SYMMETRY_CORRECTION");
    }

    // Stacking energies (bound)
    double DH_stack {0};
    double DS_stack {0};
    for (unsigned int nuc_i {0}; nuc_i != seq.size() - 1; nuc_i++) {
        string c1 {seq[nuc_i]};
        string c2 {seq[nuc_i + 1]};
        string seq_pair {c1 + c2};
        string c3 {comp_seq[nuc_i]};
        string c4 {comp_seq[nuc_i + 1]};
        string comp_seq_pair {c3 + c4};
        string key {seq_pair + "/" + comp_seq_pair};
        DH_stack += NN_Enthalpy.at(key);
        DS_stack += NN_Entropy.at(key);
    }

    // Terminal AT penalties
    int terminal_at_pairs {0};
    if (seq.front() == 'A' or seq.front() == 'T') {
        terminal_at_pairs += 1;
    }
    if (seq.back() == 'A' or seq.back() == 'T') {
        terminal_at_pairs += 1;
    }
    double DH_at {NN_Enthalpy.at("TERMINAL_AT_PENALTY") * terminal_at_pairs};
    double DS_at {NN_Entropy.at("TERMINAL_AT_PENALTY") * terminal_at_pairs};

    // Sum and correct for salt
    double DH_hybrid = DH_init + DH_stack + DH_at;
    double DS_hybrid = DS_init + DS_sym + DS_stack + DS_at;

    // Consider specifying num phosphates to account for sequences with terminal
    // residues
    DS_hybrid += 0.368 * (seq.size() / 2) * log(cation_M) / 1000;

    ThermoOfHybrid DH_DS {DH_hybrid, DS_hybrid};
    return DH_DS;
}

vector<string> find_longest_contig_complement(string seq_i, string seq_j) {
    // Find smallest sequence
    string seq_three;
    string seq_five;
    if (seq_i.size() <= seq_j.size()) {
        seq_three = seq_j;
        reverse(seq_three.begin(), seq_three.end());
        seq_five = seq_i;
    }
    else {
        seq_three = seq_i;
        reverse(seq_three.begin(), seq_three.end());
        seq_five = seq_j;
    }

    seq_three = calc_comp_seq(seq_three);

    // Iterate through all lengths and starting points
    vector<string> comp_seqs {};
    for (int subseq_len {(int)seq_three.size()}; subseq_len != 0;
         subseq_len--) {
        for (unsigned int start_i {0};
             start_i != (seq_three.size() - subseq_len + 1);
             start_i++) {
            string subseq {seq_three.substr(start_i, subseq_len)};
            size_t subseq_i {seq_five.find(subseq)};
            int subseq_count {0};
            while (subseq_i != string::npos) {
                subseq_i = seq_five.find(subseq, subseq_i + 1);
                subseq_count++;
            }
            for (int i {0}; i != subseq_count; i++) {
                comp_seqs.push_back(subseq);
            }
        }
        if (comp_seqs.empty()) {
            continue;
        }
        else {
            return comp_seqs;
        }
    }
    return comp_seqs;
}

string calc_comp_seq(string seq) {
    // Return the complementary DNA sequence.
    string comp_seq {};
    for (auto base: seq) {
        comp_seq += Complementary_Base_Pairs.at(base);
    }
    return comp_seq;
}

bool seq_is_palindromic(string seq) {
    string comp_seq {calc_comp_seq(seq)};
    string reverse_comp_seq {comp_seq};
    reverse(reverse_comp_seq.begin(), reverse_comp_seq.end());
    bool palindromic {reverse_comp_seq == seq};
    return palindromic;
}
} // namespace nearestNeighbour
