// hash.cpp

#include "hash.h"

namespace std {
    size_t hash_value(utility::VectorThree const& v) {
        size_t seed = 0;
        for (size_t i {0}; i != 3; i++) {
            hash_combine(seed, v.at(i));
        }
        return seed;
    }
}
