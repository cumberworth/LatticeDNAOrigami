// utility.cpp

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "utility.h"

namespace utility {

using std::cout;
using std::unique_ptr;
using std::vector;

int index(vector<int> container, int element) {
    for (size_t i {0}; i != container.size(); i++) {
        if (container[i] == element) {
            return i;
        }
        else
            continue;
    }
    throw NoElement {};
}

VectorThree VectorThree::operator-() {
    VectorThree neg {};
    for (size_t i {0}; i != 3; i++) {
        neg[i] = -m_container[i];
    }
    return neg;
}

VectorThree VectorThree::operator+(const VectorThree& v_2) const {
    VectorThree sum;
    for (unsigned int i {0}; i != 3; i++) {
        sum[i] = m_container[i] + v_2.at(i);
    }
    return sum;
}

VectorThree VectorThree::operator-(const VectorThree& v_2) const {
    VectorThree diff;
    for (unsigned int i {0}; i != 3; i++) {
        diff[i] = m_container[i] - v_2.at(i);
    }
    return diff;
}

bool VectorThree::operator!=(const VectorThree& v_2) const {
    for (unsigned int i {0}; i != 3; i++) {
        if (m_container[i] != v_2.at(i)) {
            return true;
        }
    }
    return false;
}

bool operator==(const VectorThree& v1, const VectorThree& v2) {
    for (unsigned int i {0}; i != 3; i++) {
        if (v1.at(i) != v2.at(i)) {
            return false;
        }
    }
    return true;
}

VectorThree VectorThree::rotate_half(VectorThree axis) {
    VectorThree rot = *this;
    axis = axis.absolute();
    if (axis == xhat) {
        rot[1] = -rot[1];
        rot[2] = -rot[2];
    }
    else if (axis == yhat) {
        rot[0] = -rot[0];
        rot[2] = -rot[2];
    }
    else if (axis == zhat) {
        rot[0] = -rot[0];
        rot[1] = -rot[1];
    }

    return rot;
}

VectorThree VectorThree::rotate(
        VectorThree origin,
        VectorThree axis,
        int turns) {

    VectorThree rot = *this;
    if (turns == 0) {
        return rot;
    }
    rot = rot - origin;
    VectorThree new_rot {rot.rotate(axis, turns)};

    return new_rot + origin;
}

VectorThree VectorThree::rotate(VectorThree axis, int turns) {
    VectorThree rot {*this};
    if (turns % 2 == 0) {
        rot = rot.rotate_half(axis);
    }
    else {
        int dir;
        bool odd_turns_even {(turns - 1) / 2 % 2 == 0};
        bool turns_neg {turns < 0};
        if (not turns_neg and odd_turns_even) {
            dir = 1;
        }
        else {
            dir = -1;
        }
        auto abs_axis {axis.absolute()};
        if (axis != abs_axis) {
            dir *= -1;
        }
        if (abs_axis == xhat) {
            std::swap(rot[1], rot[2]);
            rot[1] *= -dir;
            rot[2] *= dir;
        }
        else if (abs_axis == yhat) {
            std::swap(rot[0], rot[2]);
            rot[0] *= -dir;
            rot[2] *= dir;
        }
        else if (abs_axis == zhat) {
            std::swap(rot[0], rot[1]);
            rot[0] *= -dir;
            rot[1] *= dir;
        }
        //            else {
        //                throw SimulationMisuse {};
        //            }
    }

    return rot;
}

int VectorThree::sum() { return this->at(0) + this->at(1) + this->at(2); }

int VectorThree::abssum() {
    return abs(m_container[0]) + abs(m_container[1]) + abs(m_container[2]);
}

VectorThree VectorThree::absolute() {
    return {abs(m_container[0]), abs(m_container[1]), abs(m_container[2])};
}

VectorThree VectorThree::sort() {
    array<int, 3> sorted_container {m_container};
    std::sort(sorted_container.begin(), sorted_container.end());
    int x {sorted_container[0]};
    int y {sorted_container[1]};
    int z {sorted_container[2]};
    return {VectorThree {x, y, z}};
}

bool operator==(
        const StapleExchangeTracking& t1,
        const StapleExchangeTracking& t2) {

    return (t1.staple_insertion == t2.staple_insertion and
            t1.no_staples == t2.no_staples and
            t1.staple_type == t2.staple_type);
}

bool operator==(
        const StapleRegrowthTracking& t1,
        const StapleRegrowthTracking& t2) {

    return (t1.no_staples == t2.no_staples and
            t1.staple_type == t2.staple_type);
}

bool operator==(
        const CTCBScaffoldRegrowthTracking& t1,
        const CTCBScaffoldRegrowthTracking& t2) {

    return (t1.num_scaffold_domains == t2.num_scaffold_domains and
            t1.num_staples == t2.num_staples);
}

bool operator==(
        const ScaffoldRGRegrowthTracking& t1,
        const ScaffoldRGRegrowthTracking& t2) {

    return t1.num_scaffold_domains == t2.num_scaffold_domains;
}

bool operator==(
        const CTCBLinkerRegrowthTracking& t1,
        const CTCBLinkerRegrowthTracking& t2) {

    return (t1.central_domains_connected == t2.central_domains_connected and
            t1.num_linker_domains == t2.num_linker_domains and
            t1.num_linker_staples == t2.num_linker_staples and
            t1.num_central_domains == t2.num_central_domains and
            t1.num_central_staples == t2.num_central_staples and
            t1.disp_sum == t2.disp_sum and t1.rot_turns == t2.rot_turns);
}

// Couldn't easily template cause function changes (stoi vs stod)
vector<int> string_to_int_vector(string string_v) {
    string delim {" "};
    size_t start_pos {0};
    size_t delim_pos {string_v.find(delim)};
    vector<int> v {};
    while (delim_pos != string::npos) {
        v.push_back(stoi(string_v.substr(start_pos, delim_pos)));
        start_pos = delim_pos + 1;
        delim_pos = string_v.find(delim, start_pos);
    }
    v.push_back(stod(string_v.substr(start_pos, string_v.size())));

    return v;
}

vector<double> string_to_double_vector(string string_v) {
    string delim {" "};
    size_t start_pos {0};
    size_t delim_pos {string_v.find(delim)};
    vector<double> v {};
    while (delim_pos != string::npos) {
        v.push_back(stod(string_v.substr(start_pos, delim_pos)));
        start_pos = delim_pos + 1;
        delim_pos = string_v.find(delim, start_pos);
    }
    v.push_back(stod(string_v.substr(start_pos, string_v.size())));

    return v;
}

vector<string> string_to_string_vector(string string_v) {
    std::stringstream sstream {string_v};
    vector<string> v {};
    while (not sstream.eof()) {
        string ele;
        sstream >> ele;
        v.push_back(ele);
    }

    return v;
}

template <typename Out>
void split(const string& s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

vector<string> split(const string& s, char delim) {
    vector<string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

Fraction::Fraction(string unparsed_fraction) {
    string delimiter {"/"};
    auto delim_pos {unparsed_fraction.find(delimiter)};

    // Assume is real number
    if (delim_pos == string::npos) {
        m_numerator = stod(unparsed_fraction);
        m_denominator = 1;
    }
    else {
        auto end_pos {unparsed_fraction.size()};
        m_numerator = stod(unparsed_fraction.substr(0, delim_pos));
        m_denominator = stod(unparsed_fraction.substr(delim_pos + 1, end_pos));
    }
    m_double_fraction = m_numerator / m_denominator;
}

} // namespace utility
