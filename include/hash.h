// hash.h
#ifndef HASH_H
#define HASH_H

#include <iostream>

#include "utility.h"

namespace std {

    using utility::VectorThree;
    using utility::StapleExchangeTracking;
    using utility::StapleRegrowthTracking;
    using utility::CTCBScaffoldRegrowthTracking;
    using utility::CTCBLinkerRegrowthTracking;

    /* Copied from a stack exchange question (which is copied from the BOOST
       library) for allowing pairs to be hashed.
    */
    template <class T>
    inline void hash_combine(std::size_t & seed, const T & v)
    {
      hash<T> hasher;
      seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    template<typename S, typename T> struct hash<pair<S, T>> {
        inline size_t operator()(const pair<S, T>& v) const {
            size_t seed = 0;
            hash_combine(seed, v.first);
            hash_combine(seed, v.second);
            return seed;
        }
    };

    template<typename T> struct hash<vector<T>> {
        inline size_t operator()(const vector<T>& v) const {
            size_t seed = 0;
            for (auto i: v) {
                hash_combine(seed, i);
            }
            return seed;
        }
    };

    template<> struct hash<VectorThree> {
        inline size_t operator()(const VectorThree& v) const {
            size_t seed = 0;
            for (size_t i {0}; i != 3; i++) {
                hash_combine(seed, v.at(i));
            }
            return seed;
        }
    };

    template<> struct hash<StapleExchangeTracking> {
        inline size_t operator()(const StapleExchangeTracking& x) const {
            size_t seed = 0;
            hash_combine(seed, x.no_staples);
            hash_combine(seed, x.staple_type);
            return seed;
        }
    };

    template<> struct hash<StapleRegrowthTracking> {
        inline size_t operator()(const StapleRegrowthTracking& x) const {
            size_t seed = 0;
            hash_combine(seed, x.no_staples);
            hash_combine(seed, x.staple_type);
            return seed;
        }
    };

    template<> struct hash<CTCBScaffoldRegrowthTracking> {
        inline size_t operator()(const CTCBScaffoldRegrowthTracking& x) const {
            size_t seed = 0;
            hash_combine(seed, x.num_scaffold_domains);
            hash_combine(seed, x.num_staples);
            return seed;
        }
    };

    template<> struct hash<CTCBLinkerRegrowthTracking> {
        inline size_t operator()(const CTCBLinkerRegrowthTracking& x) const {
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

}

#endif // HASH_H
