// hash.h
#ifndef HASH_H
#define HASH_H

#include <iostream>

#include "boost/functional/hash.hpp"

#include "utility.h"

namespace std {

using boost::hash_combine;
using boost::hash_value;

template <>
struct hash<utility::VectorThree> {
    inline size_t operator()(const utility::VectorThree& v) const {
        size_t seed = 0;
        for (size_t i {0}; i != 3; i++) {
            hash_combine(seed, v.at(i));
        }
        return seed;
    }
};

size_t hash_value(utility::VectorThree const& v);

template <>
struct hash<utility::StapleExchangeTracking> {
    inline size_t operator()(const utility::StapleExchangeTracking& x) const {
        size_t seed = 0;
        hash_combine(seed, x.no_staples);
        hash_combine(seed, x.staple_type);
        return seed;
    }
};

template <>
struct hash<utility::StapleRegrowthTracking> {
    inline size_t operator()(const utility::StapleRegrowthTracking& x) const {
        size_t seed = 0;
        hash_combine(seed, x.no_staples);
        hash_combine(seed, x.staple_type);
        return seed;
    }
};

template <>
struct hash<utility::ScaffoldRGRegrowthTracking> {
    inline size_t operator()(
            const utility::ScaffoldRGRegrowthTracking& x) const {
        size_t seed = 0;
        hash_combine(seed, x.num_scaffold_domains);
        return seed;
    }
};

template <>
struct hash<utility::CTCBScaffoldRegrowthTracking> {
    inline size_t operator()(
            const utility::CTCBScaffoldRegrowthTracking& x) const {
        size_t seed = 0;
        hash_combine(seed, x.num_scaffold_domains);
        hash_combine(seed, x.num_staples);
        return seed;
    }
};

template <>
struct hash<utility::CTCBLinkerRegrowthTracking> {
    inline size_t operator()(
            const utility::CTCBLinkerRegrowthTracking& x) const {
        size_t seed = 0;
        hash_combine(seed, x.central_domains_connected);
        hash_combine(seed, x.num_linker_domains);
        hash_combine(seed, x.num_linker_staples);
        hash_combine(seed, x.num_central_domains);
        hash_combine(seed, x.num_central_staples);
        hash_combine(seed, x.disp_sum);
        hash_combine(seed, x.rot_turns);
        return seed;
    }
};

template <typename S, typename T>
struct hash<pair<S, T>> {
    inline size_t operator()(const pair<S, T>& v) const {
        size_t seed {0};
        boost::hash_combine(seed, hash_value(v.first));
        boost::hash_combine(seed, hash_value(v.second));
        return seed;
    }
};

template <typename T>
struct hash<vector<T>> {
    inline size_t operator()(const vector<T>& v) const {
        size_t seed = 0;
        for (auto i: v) {
            hash_combine(seed, i);
        }
        return seed;
    }
};

} // namespace std

#endif // HASH_H
