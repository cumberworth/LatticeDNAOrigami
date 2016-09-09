// ideal_random_walk.cpp

#include <algorithm>
#include <boost/math/special_functions/factorials.hpp>

#include "ideal_random_walk.h"

using boost::math::factorial;

using namespace IdealRandomWalk;

unsigned long long IdealRandomWalks::num_walks(
        VectorThree start_pos,
        VectorThree end_pos,
        int steps) {

    // Check stored values
    VectorThree DR {end_pos - start_pos};
    pair<VectorThree, int> walk_key {DR, steps};
    try {
        return m_num_walks.at(walk_key);
    }
    catch (std::out_of_range) {
    }
    int DX {DR[0]};
    int DY {DR[0]};
    int DZ {DR[0]};
    int Nminus {(steps - DX - DY - DZ) / 2};
    int Nplus {(steps - DX - DY + DZ) / 2};
    unsigned long long walks;
    
    // Need some negative steps to reach a negative location
	for (int ybar {0}; ybar != Nminus + 1; ybar++) {
        for (int xbar {0}; xbar != Nminus + 1 - ybar; xbar++) {
            auto f1 {factorial<unsigned long long>(steps)};
            auto f2 {factorial<unsigned long long>(xbar)};
            if (xbar + DX < 0) {
                continue;
            }

            auto f3 {factorial<unsigned long long>(xbar + DX)};
            auto f4 {factorial<unsigned long long>(ybar)};
            if (ybar + DY < 0) {
                continue;
            }

            auto f5 {factorial<unsigned long long>(ybar + DY)};
            auto f6 {factorial<unsigned long long>(Nminus - xbar - ybar)};
            if (Nplus - xbar - ybar < 0) {
                continue;
            }

            auto f7 {factorial<unsigned long long>(Nplus - xbar - ybar)};
            walks += f1 / (f2 * f3 * f4 * f5 * f6 * f7);
        }
    }

    // Add entries
    m_num_walks[walk_key] =  walks;

    return walks;
}
