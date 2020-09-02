// calc_num_walks.h

#ifndef CALC_NUM_WALKS_H
#define CALC_NUM_WALKS_H

#include <unordered_map>
#include <utility>
#include <string>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/array.hpp>

#include "parser.h"
#include "ideal_random_walk.h"

using namespace idealRandomWalk;
using namespace parser;

namespace calcNumWalks {

    void calc_num_ideal_walks(int max_d, int max_N, string filename);
}

#endif // CALC_NUM_WALKS_H
