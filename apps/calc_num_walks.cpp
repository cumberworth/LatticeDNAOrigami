// calc_num_walks.cpp

#include <fstream>
#include <string>
#include <unordered_map>
#include <utility>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/utility.hpp>

#include "LatticeDNAOrigami/ideal_random_walk.hpp"
#include "LatticeDNAOrigami/parser.hpp"
#include "LatticeDNAOrigami/utility.hpp"

using namespace idealRandomWalk;

int main(int argc, char* argv[]) {
    int max_d {78};
    int max_N {156};
    string filename {"num_walks.arch"};
    IdealRandomWalks ideal_random_walks {};
    VectorThree start_pos {0, 0, 0};
    for (int x {0}; x <= max_d; x++) {
        std::cout << x << "\n";
        for (int y {0}; y <= x; y++) {
            for (int z {0}; z <= y; z++) {
                VectorThree end_pos {x, y, z};
                for (int N {0}; N <= max_N; N++) {
                    auto walks {ideal_random_walks.num_walks(
                            start_pos, end_pos, N)};
                    if (walks == 0) {
                        ideal_random_walks.delete_entry(start_pos, end_pos, N);
                    }
                }
            }
        }
    }
    std::ofstream file {filename};
    boost::archive::binary_oarchive arch {file};
    arch << ideal_random_walks;
}
