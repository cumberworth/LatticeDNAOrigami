// ideal_random_walk.h

#ifndef IDEAL_RANDOM_WALK_H
#define IDEAL_RANDOM_WALK_H

#include <unordered_map>
#include <utility>

#include "hash.h"

using std::unordered_map;
using std::pair;

namespace IdealRandomWalk {

    class IdealRandomWalks {
        public:
            unsigned long long num_walks(VectorThree start_pos,
                    VectorThree end_pos, int steps);
            
        private:
            unordered_map<pair<VectorThree, int>, unsigned long long> m_num_walks {};
    };
}

#endif // IDEAL_RANDOM_WALK_H
