// ideal_random_walk.cpp

#include <algorithm>
#include <boost/math/special_functions/factorials.hpp>

#include "ideal_random_walk.h"

using boost::math::factorial;

using namespace IdealRandomWalk;

double IdealRandomWalks::num_walks(
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
    int DY {DR[1]};
    int DZ {DR[2]};
    //cout << "Steps: " << steps << "\n";
    //cout << "DR: " << DX << " " << DY << " " << DZ << "\n";
    int Nminus {(steps - DX - DY - DZ) / 2};
    //cout << "Nminus: " << Nminus << "\n";
    int Nplus {(steps - DX - DY + DZ) / 2};
    double walks {0};

    int DR_sum {DX + DY + DZ};
    if (DR_sum > steps or (steps - DR_sum) % 2 != 0) {

        // Add entry
        m_num_walks[walk_key] =  walks;
        //cout << "Number of walks: " << walks << "\n\n";
        return walks;
    }
    
    // Need some negative steps to reach a negative location
	for (int ybar {0}; ybar != Nminus + 1; ybar++) {
        for (int xbar {0}; xbar != Nminus + 1 - ybar; xbar++) {
            auto f1 {factorial<double>(steps)};
            auto f2 {factorial<double>(xbar)};
            if (xbar + DX < 0) {
                continue;
            }

            auto f3 {factorial<double>(xbar + DX)};
            auto f4 {factorial<double>(ybar)};
            if (ybar + DY < 0) {
                continue;
            }

            auto f5 {factorial<double>(ybar + DY)};
            auto f6 {factorial<double>(Nminus - xbar - ybar)};
            if (Nplus - xbar - ybar < 0) {
                continue;
            }

            auto f7 {factorial<double>(Nplus - xbar - ybar)};
            walks += f1 / (f2 * f3 * f4 * f5 * f6 * f7);
        }
    }

    // Add entries
    m_num_walks[walk_key] =  walks;
    //cout << "Number of walks: " << walks << "\n\n";

    return walks;
}
