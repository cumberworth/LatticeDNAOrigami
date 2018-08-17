// hash.h
#ifndef HASH_H
#define HASH_H

#include <iostream>

#include "boost/functional/hash.hpp"

#include "utility.h"

namespace std {

    template<typename S, typename T> struct hash<pair<S, T>> {
        inline size_t operator()(const pair<S, T>& v) const {
            size_t seed = 0;
            boost::hash_combine(seed, v.first);
            boost::hash_combine(seed, v.second);
            return seed;
        }
    };

    template<typename T> struct hash<vector<T>> {
        inline size_t operator()(const vector<T>& v) const {
            size_t seed = 0;
            for (auto i: v) {
                boost::hash_combine(seed, i);
            }
            return seed;
        }
    };

    template<> struct hash<utility::VectorThree> {
        inline size_t operator()(const utility::VectorThree& v) const {
            size_t seed = 0;
            for (size_t i {0}; i != 3; i++) {
                boost::hash_combine(seed, v.at(i));
            }
            return seed;
        }
    };

    template<> struct hash<utility::StapleExchangeTracking> {
        inline size_t operator()(const utility::StapleExchangeTracking& x) const {
            size_t seed = 0;
            boost::hash_combine(seed, x.no_staples);
            boost::hash_combine(seed, x.staple_type);
            return seed;
        }
    };

    template<> struct hash<utility::StapleRegrowthTracking> {
        inline size_t operator()(const utility::StapleRegrowthTracking& x) const {
            size_t seed = 0;
            boost::hash_combine(seed, x.no_staples);
            boost::hash_combine(seed, x.staple_type);
            return seed;
        }
    };

    template<> struct hash<utility::ScaffoldRGRegrowthTracking> {
        inline size_t operator()(const utility::ScaffoldRGRegrowthTracking& x) const {
            size_t seed = 0;
            boost::hash_combine(seed, x.num_scaffold_domains);
            return seed;
        }
    };

    template<> struct hash<utility::CTCBScaffoldRegrowthTracking> {
        inline size_t operator()(const utility::CTCBScaffoldRegrowthTracking& x) const {
            size_t seed = 0;
            boost::hash_combine(seed, x.num_scaffold_domains);
            boost::hash_combine(seed, x.num_staples);
            return seed;
        }
    };

    template<> struct hash<utility::CTCBLinkerRegrowthTracking> {
        inline size_t operator()(const utility::CTCBLinkerRegrowthTracking& x) const {
            size_t seed = 0;
            boost::hash_combine(seed, x.central_domains_connected);
            boost::hash_combine(seed, x.num_linker_domains);
            boost::hash_combine(seed, x.num_linker_staples);
            boost::hash_combine(seed, x.num_central_domains);
            boost::hash_combine(seed, x.num_central_staples);
            boost::hash_combine(seed, x.disp_sum);
            boost::hash_combine(seed, x.rot_turns);
            return seed;
        }
    };

}

#endif // HASH_H
