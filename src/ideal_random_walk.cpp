// ideal_random_walk.cpp

#include <cstdint>

#include <boost/math/special_functions/factorials.hpp>

#include "LatticeDNAOrigami/ideal_random_walk.hpp"

namespace idealRandomWalk {

using boost::math::factorial;

long double IdealRandomWalks::num_walks(
        VectorThree start_pos,
        VectorThree end_pos,
        int steps) {

    // The count only depends on one permutation of the absolute displacement
    VectorThree DR {(end_pos - start_pos).absolute().sort()};
    if (steps < 0 or steps > 0xFFFF or DR[2] > 0xFFFF) {
        throw utility::OrigamiMisuse {"Random walk outside supported range"};
    }
    std::uint64_t key {
            static_cast<std::uint64_t>(DR[0]) << 48 |
            static_cast<std::uint64_t>(DR[1]) << 32 |
            static_cast<std::uint64_t>(DR[2]) << 16 |
            static_cast<std::uint64_t>(steps)};
    auto stored {m_num_walks.find(key)};
    if (stored != m_num_walks.end()) {
        return stored->second;
    }
    long double walks {calc_num_walks(DR[0], DR[1], DR[2], steps)};
    m_num_walks.emplace(key, walks);

    return walks;
}

long double IdealRandomWalks::calc_num_walks(
        int DX,
        int DY,
        int DZ,
        int steps) {
    int Nminus {(steps - DX - DY - DZ) / 2};
    int Nplus {(steps - DX - DY + DZ) / 2};
    long double walks {0};

    int DR_sum {DX + DY + DZ};
    if (DR_sum > steps or (steps - DR_sum) % 2 != 0) {
        return walks;
    }

    // Need some negative steps to reach a negative location
    for (int ybar {0}; ybar != Nminus + 1; ybar++) {
        for (int xbar {0}; xbar != Nminus + 1 - ybar; xbar++) {
            auto f1 {factorial<long double>(steps)};
            auto f2 {factorial<long double>(xbar)};
            if (xbar + DX < 0) {
                continue;
            }

            auto f3 {factorial<long double>(xbar + DX)};
            auto f4 {factorial<long double>(ybar)};
            if (ybar + DY < 0) {
                continue;
            }

            auto f5 {factorial<long double>(ybar + DY)};
            auto f6 {factorial<long double>(Nminus - xbar - ybar)};
            if (Nplus - xbar - ybar < 0) {
                continue;
            }

            auto f7 {factorial<long double>(Nplus - xbar - ybar)};
            walks += f1 / (f2 * f3 * f4 * f5 * f6 * f7);
        }
    }

    return walks;
}
} // namespace idealRandomWalk
