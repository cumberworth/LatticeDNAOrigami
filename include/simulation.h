// simulation.h

#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <memory>
#include <map>
#include <iostream>
#include <utility>
#include <set>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/process.hpp>

#include "parser.h"
#include "origami_system.h"
#include "movetypes.h"
#include "files.h"
#include "ideal_random_walk.h"
#include "order_params.h"

namespace Simulation {

    using std::vector;
    using std::unique_ptr;
    using std::shared_ptr;
    using std::map;
    using std::ostream;
    using std::ofstream;
    using std::unordered_map;
    using std::set;

    namespace mpi = boost::mpi;
    namespace bp = boost::process;

    using namespace Parser;
    using namespace Origami;
    using namespace Movetypes;
    using namespace Files;
    using namespace IdealRandomWalk;
    using namespace OrderParams;

    vector<OrigamiOutputFile*> setup_output_files(
            InputParameters& params,
            string output_filebase,
            OrigamiSystem& origami);

    class GCMCSimulation {
        public:
            GCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            virtual ~GCMCSimulation();
            virtual void run() = 0;

        protected:
            OrigamiSystem& m_origami_system;
            ostream* m_logging_stream;
            int m_logging_freq;
            int m_centering_freq;
            int m_constraint_check_freq;
            int m_vmd_pipe_freq;
            InputParameters& m_params;
            vector <OrigamiOutputFile*> m_output_files;
            vector<shared_ptr<MCMovetype>> m_movetypes;
            vector<double> m_cumulative_probs;
            RandomGens m_random_gens {};
            IdealRandomWalks m_ideal_random_walks {};

            // VMD realtime visualization
            OrigamiVSFOutputFile* vmd_struct_file {NULL};
            OrigamiVCFOutputFile* vmd_coors_file {NULL};
            OrigamiStateOutputFile* vmd_states_file {NULL};
            OrigamiOrientationOutputFile* vmd_ores_file {NULL};
            bp::child* vmd_proc {NULL};

            // Shared interface
            virtual void update_internal(long long int step) = 0;

            // Shared methods
            void construct_movetypes(InputParameters& params);
            void simulate(long long int steps, long long int start_step=0);
            shared_ptr<MCMovetype> select_movetype();
            void write_log_entry(
                    long long int step,
                    MCMovetype& movetype,
                    bool accepted);
            void setup_vmd_pipe();
            void pipe_to_vmd();
            void close_vmd_pipe();
            void close_output_files();
    };
}

#endif // SIMULATION_H
