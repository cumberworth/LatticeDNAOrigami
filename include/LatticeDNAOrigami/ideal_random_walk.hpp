// ideal_random_walk.h

#ifndef IDEAL_RANDOM_WALK_H
#define IDEAL_RANDOM_WALK_H

#include <string>
#include <unordered_map>
#include <utility>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/utility.hpp>

#include "LatticeDNAOrigami/hash.hpp"
#include "LatticeDNAOrigami/utility.hpp"

namespace idealRandomWalk {

using std::pair;
using std::string;
using std::unordered_map;

using utility::VectorThree;

class IdealRandomWalks {
  public:
    long double num_walks(
            VectorThree start_pos,
            VectorThree end_pos,
            int steps);

    void delete_entry(VectorThree start_pos, VectorThree end_pos, int steps);

  private:
    unordered_map<pair<VectorThree, int>, long double> m_num_walks {};
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& arch, const unsigned int) {
        arch& m_num_walks;
    }
};
} // namespace idealRandomWalk

#endif // IDEAL_RANDOM_WALK_H
