// random_gens.h

#ifndef RANDOM_GENS_H
#define RANDOM_GENS_H

#include <random>
#include <unordered_map>
#include <utility>

#include "hash.h"

namespace randomGen {

using std::pair;
using std::unordered_map;

class RandomGens {
  public:
    RandomGens();
    ~RandomGens();

    void set_seed(int seed);
    double uniform_real();
    int uniform_int(int lower, int upper);

  private:
    std::mt19937_64 m_random_engine {};
    std::uniform_real_distribution<double> m_uniform_real_dist;
    unordered_map<pair<int, int>, std::uniform_int_distribution<int>&>
            m_uniform_int_dists {};
};
} // namespace randomGen

#endif // RANDOM_GENS_H
