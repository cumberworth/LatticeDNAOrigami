// main.cpp

#include <iostream>

#include "parser.h"
#include "nearest_neighbour.h"
#include "domain.h"
#include "origami_system.h"
#include "files.h"

using std::cout;

using namespace Utility;
using namespace Files;
using namespace Origami;
using namespace DomainContainer;

int main() {
    OrigamiInputFile origami_input {"tests/snodin_assembled.json"};
    vector<vector<int>> identities {origami_input.m_identities};
    vector<vector<string>> sequences {origami_input.m_sequences};
    vector<Chain> chains {origami_input.m_chains};
    OrigamiSystem origami {
            identities,
            sequences,
            chains,
            300,
            1,
            1,
            1,
            false};
    //origami.add_chain(1);
    //Domain& cd_i {*origami.m_domains[1][0]};
    //VectorThree pos_1 {0, 0, 0};
    //VectorThree ore_1 {-1, 0, 0};
    //origami.set_domain_config(cd_i, pos_1, ore_1);
    //VectorThree pos_2 {1, 0, 0};
    //VectorThree ore_2 {1, 0, 0};
    //origami.set_domain_config(cd_i, pos_2, ore_2);
    // parse input file
    // setup io
    // setup sim
    // run sim
}
