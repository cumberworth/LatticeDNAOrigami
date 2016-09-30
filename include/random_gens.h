// random_gens.h

#ifndef RANDOM_GENS_H
#define RANDOM_GENS_H

#include <random>
#include <utility>
#include <unordered_map>

#include "hash.h"

using std::pair;
using std::unordered_map;

namespace RandomGen {
    class RandomGens {
        public:
            RandomGens();
            ~RandomGens();

            std::mt19937_64 m_random_engine {};
            std::uniform_real_distribution<double> m_uniform_real_dist;
            unordered_map<pair<int, int>, std::uniform_int_distribution<int>&> m_uniform_int_dists {};

            double uniform_real() {return m_uniform_real_dist(m_random_engine);};
            int uniform_int(int lower, int upper);
    };
}

#endif // RANDOM_GENS_H
