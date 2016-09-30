// random_gens.cpp

#include <random>
#include <iostream>

#include "random_gens.h"

using namespace RandomGen;
using std::cout;

RandomGens::RandomGens() {

    // Seed random number generator 
    std::random_device true_random_engine {}; 
    auto seed {true_random_engine()}; 
    m_random_engine.seed(seed); 
}

RandomGens::~RandomGens() {
    for (auto const key: m_uniform_int_dists) {
        delete &m_uniform_int_dists.at(key.first);
    }
}

int RandomGens::uniform_int(int lower, int upper) {

    // Check if distribution used previously
    pair<int, int> key {lower, upper};
    try {
        std::uniform_int_distribution<int>& dist {m_uniform_int_dists.at(key)};
        return dist(m_random_engine);
    }

    // If not, make it and store it
    catch (std::out_of_range) {
        std::uniform_int_distribution<int>* dist {new std::uniform_int_distribution<int> {lower, upper}};
        int random_int {(*dist)(m_random_engine)};
        m_uniform_int_dists.insert({key, *dist});
        return random_int;
    }
}
