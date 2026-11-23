// ideal_random_walk.h

#ifndef IDEAL_RANDOM_WALK_H
#define IDEAL_RANDOM_WALK_H

#include <cstdint>

#include <boost/unordered/unordered_flat_map.hpp>

#include "LatticeDNAOrigami/utility.hpp"

namespace idealRandomWalk {

using utility::VectorThree;

// Number of ideal (non-self-avoiding) random walks on a cubic lattice between
// two points with a given number of steps, computed on demand and memoized
class IdealRandomWalks {
  public:
    long double num_walks(
            VectorThree start_pos,
            VectorThree end_pos,
            int steps);

  private:
    static long double calc_num_walks(int dx, int dy, int dz, int steps);

    // Keyed on the sorted absolute displacement and the number of steps,
    // packed into 16 bits each
    boost::unordered_flat_map<std::uint64_t, long double> m_num_walks {};
};
} // namespace idealRandomWalk

#endif // IDEAL_RANDOM_WALK_H
