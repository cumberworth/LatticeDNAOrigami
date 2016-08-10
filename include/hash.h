// hash.h
#ifndef HASH_H
#define HASH_H

#include "utility.h"

using Utility::CDPair;
using Utility::VectorThree;

/* Copied from a stack exchange question (which is copied from the BOOST
   library) for allowing pairs to be hashed.
*/
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {

    template<typename S, typename T> struct hash<pair<S, T>> {
        inline size_t operator()(const pair<S, T>& v) const {
            size_t seed = 0;
            ::hash_combine(seed, v.first);
            ::hash_combine(seed, v.second);
            return seed;
        }
    };

    template<typename T> struct hash<vector<T>> {
        inline size_t operator()(const vector<T>& v) const {
            size_t seed = 0;
            for (auto i: v) {
                ::hash_combine(seed, i);
            }
            return seed;
        }
    };

    template<> struct hash<VectorThree> {
        inline size_t operator()(const VectorThree& v) const {
            size_t seed = 0;
            for (unsigned int i {0}; i != 3; i++) {
                ::hash_combine(seed, v.at(i));
            }
            return seed;
        }
    };

    template<> struct hash<CDPair> {
        inline size_t operator()(const CDPair& v) const {
            size_t seed = 0;
            ::hash_combine(seed, v.c);
            ::hash_combine(seed, v.d);
            return seed;
        }
    };
}

#endif // HASH_H
