// calc_num_walks.cpp

#include <fstream>

#include "calc_num_walks.h"
#include "parser.h"
#include "utility.h"

using namespace CalcNumWalks;
using namespace Parser;

int main(int argc, char* argv[]) {
    InputParameters params {argc, argv};
    calc_num_ideal_walks(78, 156, params.m_num_walks_filename);
}

void CalcNumWalks::calc_num_ideal_walks(int max_d, int max_N, string filename) {
    IdealRandomWalks ideal_random_walks {};
    VectorThree start_pos {0, 0, 0};
    for (int x {0}; x <= max_d; x++) {
        cout << x << "\n";
        for (int y {0}; y <= x; y++) {
            for (int z {0}; z <= y; z++) {
                VectorThree end_pos {x, y, z};
                for (int N {0}; N <= max_N; N++) {
                    auto walks {ideal_random_walks.num_walks(start_pos, end_pos, N)};
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
