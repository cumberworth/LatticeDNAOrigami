// random_gens.cpp

#include <iostream>
#include <random>

#include "LatticeDNAOrigami/random_gens.hpp"

namespace randomGen {

using std::cout;

RandomGens::RandomGens() {

    // Seed random number generator
    std::random_device true_random_engine {};
    auto seed {true_random_engine()};
    cout << "Truly random seed: " << seed << "\n";
    m_random_engine.seed(seed);
}

void RandomGens::set_seed(int seed) { m_random_engine.seed(seed); }

double RandomGens::uniform_real() {
    return m_uniform_real_dist(m_random_engine);
}

int RandomGens::uniform_int(int lower, int upper) {

    // Constructing an integer distribution is cheap and it holds no state
    // between calls, so there is no need to cache one per range
    return std::uniform_int_distribution<int> {lower, upper}(m_random_engine);
}
} // namespace randomGen
