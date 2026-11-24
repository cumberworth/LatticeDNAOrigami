// hash.h
#ifndef HASH_H
#define HASH_H

#include <cstdint>
#include <iostream>

#include <boost/unordered/unordered_flat_map.hpp>

#include "boost/functional/hash.hpp"

#include "LatticeDNAOrigami/utility.hpp"

namespace std {

using boost::hash_combine;
using boost::hash_value;

// Multiply each coordinate by a different odd constant and combine. The flat
// hash maps that use this mix the result further, so it need not avalanche.
template <>
struct hash<utility::VectorThree> {
    inline size_t operator()(const utility::VectorThree& v) const {
        auto x {static_cast<uint64_t>(static_cast<uint32_t>(v[0]))};
        auto y {static_cast<uint64_t>(static_cast<uint32_t>(v[1]))};
        auto z {static_cast<uint64_t>(static_cast<uint32_t>(v[2]))};
        return static_cast<size_t>(
                x * 0x9E3779B97F4A7C15ull ^ y * 0xC2B2AE3D27D4EB4Full ^
                z * 0x165667B19E3779F9ull);
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

namespace hashing {

// Hash for (chain index, domain index) style keys: pack both into 64 bits.
// The flat maps mix the result further.
struct IntPairHash {
    size_t operator()(const std::pair<int, int>& p) const noexcept {
        return static_cast<size_t>(
                static_cast<uint64_t>(static_cast<uint32_t>(p.first)) << 32 |
                static_cast<uint32_t>(p.second));
    }
};

// Open-addressing map keyed on a pair of ints. Unlike std::unordered_map,
// inserting can move elements, so don't hold references across inserts.
template <typename V>
using IntPairMap =
        boost::unordered_flat_map<std::pair<int, int>, V, IntPairHash>;

} // namespace hashing

#endif // HASH_H
