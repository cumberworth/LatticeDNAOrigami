// simulation.cpp

#include <random>
#include <iostream>

#include "utility.h"
#include "random_gens.h"
#include "movetypes.h"
#include "simulation.h"

using std::cout;

using namespace Movetypes;
using namespace Simulation;
using namespace Utility;
using namespace RandomGen;

GCMCSimulation::GCMCSimulation(
        OrigamiSystem origami_system,
        OrigamiTrajOutputFile output_file,
        vector<MovetypeConstructor> movetype_constructors,
        vector<double> movetype_probs) :
        m_origami_system {origami_system}, m_output_file {output_file}, m_movetype_constructors {movetype_constructors} {

    // Create cumulative probability array
    double cum_prob {0};
    for (size_t i {0}; i != m_movetype_constructors.size(); i++) {
        cum_prob += movetype_probs[i];
        m_cumulative_probs.push_back(cum_prob);
    }
}

void GCMCSimulation::run(int steps, int logging_freq=1, int center_freq=1) {
    for (int step {1}; step != (steps + 1); step ++) {
        unique_ptr<MCMovetype> movetype {select_movetype()};
        bool accepted;
        try {
            accepted = movetype->attempt_move();
        }
        catch (MoveRejection) {
            accepted = false;
        }

        if (not accepted) {
            movetype->reset_origami();
        }

        if (step != center_freq) {
            m_origami_system.centre();
            m_origami_system.check_all_constraints();
        }

        if (step % logging_freq == 0) {
            write_log_entry(step, accepted);
        }

        // Write to traj file(s)
    }
}

unique_ptr<MCMovetype> GCMCSimulation::select_movetype() {
    unique_ptr<MCMovetype> movetype;
    double prob {m_random_gens.uniform_real()};
    for (size_t i {0}; i != m_cumulative_probs.size(); i++) {
        if (prob < m_cumulative_probs[i]) {
            movetype = m_movetype_constructors[i](m_origami_system, m_random_gens);
            break;
        }
    }
    return movetype;
}

void GCMCSimulation::write_log_entry(int step, bool accepted) {
    cout << "Step: " << step << "\n";
    cout << "Accepted: " << accepted << "\n";
    cout << "\n";
}
