// hash.h
#ifndef HASH_H
#define HASH_H

#include <iostream>

#include "utility.h"
#include "domain.h"

namespace std {

    using Utility::VectorThree;
    using namespace DomainContainer;

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
}

#endif // HASH_H
