// simulation.h

#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <memory>
#include <map>
#include <iostream>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "parser.h"
#include "origami_system.h"
#include "movetypes.h"
#include "files.h"
#include "ideal_random_walk.h"

using std::vector;
using std::unique_ptr;
using std::map;
using std::ostream;

namespace mpi = boost::mpi;

using namespace Parser;
using namespace Origami;
using namespace Movetypes;
using namespace Files;
using namespace IdealRandomWalk;

namespace Simulation {

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
            InputParameters& m_params;
            vector <OrigamiOutputFile*> m_output_files;
            vector<MovetypeConstructor> m_movetype_constructors;
            vector<double> m_cumulative_probs;
            RandomGens m_random_gens {};
            IdealRandomWalks m_ideal_random_walks {};

            // Shared methods
            void simulate(long int steps, int start_step=0);
            unique_ptr<MCMovetype> select_movetype();
            void write_log_entry(
                    int step,
                    MCMovetype& movetype,
                    bool accepted);
    };

    class ConstantTGCMCSimulation: public GCMCSimulation {
        public:
            ConstantTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            void run() {simulate(m_steps);}
        private:
            long int m_steps;
    };

    class AnnealingGCMCSimulation: public GCMCSimulation {
        public:
            AnnealingGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            void run();
        private:
            double m_max_temp;
            double m_min_temp;
            double m_temp_interval;
            int m_steps_per_temp;
    };

    class PTGCMCSimulation: public GCMCSimulation {
        public:
            PTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            ~PTGCMCSimulation() {delete m_logging_stream;}
            void run();

        protected:
            mpi::environment m_env;
            mpi::communicator m_world;
            int m_rank {m_world.rank()};
            int m_master_rep {0};
            double m_temp;
            double m_staple_u;
            double m_volume;
            int m_swaps;
            int m_num_reps;
            int m_exchange_interval;

            // Data only filled in master
            ofstream m_swapfile;
            vector<double> m_temps;
            vector<double> m_staple_us;
            vector<int> m_tempi_to_repi;

            void slave_send_and_recieve(int swap_i);
            void master_receive(
                        int swap_i,
                        vector<double>& enthalpies,
                        vector<int>& staples);
            void master_send(int swap_i);
            void attempt_exchange(
                    int swap_i,
                    vector<int>& attempt_count,
                    vector<int>& swap_count);
            bool test_acceptance(
                    double temp1,
                    double temp2,
                    double staple_u1,
                    double staple_u2,
                    double enthalpy1,
                    double enthalpy2,
                    int N1,
                    int N2);
            void write_swap_entry();
            void write_acceptance_freqs(
                    vector<int> attempt_count,
                    vector<int> swap_count);
    };
}

#endif // SIMULATION_H
