// random_gens.h

#ifndef RANDOM_GENS_H
#define RANDOM_GENS_H

#include <random>

namespace randomGen {

class RandomGens {
  public:
    RandomGens();

    void set_seed(int seed);
    double uniform_real();
    int uniform_int(int lower, int upper);
    std::mt19937_64 m_random_engine {};

  private:
    std::uniform_real_distribution<double> m_uniform_real_dist;
};
} // namespace randomGen

#endif // RANDOM_GENS_H
