#include <catch.hpp>

#define private public
#define protected public

#include "parser.h"
#include "nearest_neighbour.h"
#include "domain.h"
#include "files.h"
#include "simulation.h"

using std::cout;

using namespace Parser;
using namespace Utility;
using namespace Files;
using namespace Origami;
using namespace DomainContainer;
using namespace Simulation;
using namespace Movetypes;
using namespace NearestNeighbour;

SCENARIO("Simulation methods are run") {
    //setup shared objects
    GIVEN("Constant T simulation is run") {
        //very simple tests that it runs
    }
    GIVEN("T annealing simulation is run") {
        //check that t actually increments correctly
        //check start and end points
    }
    GIVEN("PTMC simulation is run") {
        //check that exchange probability calculated matches expected with set t, u, n, and e diff
    }
}
