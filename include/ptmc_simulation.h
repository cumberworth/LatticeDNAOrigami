// ptmc_simulation.h

#ifndef PTMC_SIMULATION_H
#define PTMC_SIMULATION_H

#include "simulation.h"

using namespace Simulation;

namespace PTMC {

    class PTGCMCSimulation: public GCMCSimulation {
        // Base method for parallel tempering in GC ensemble
        // Exchange in only 1D, but can exchange any number of quantities
        public:
            PTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            ~PTGCMCSimulation();
            void run();

        protected:

            // MPI variables
            mpi::environment m_env;
            mpi::communicator m_world;

            // General PTMC parameters
            int m_rank {m_world.rank()};
            int m_master_rep {0};
            int m_swaps;
            int m_num_reps;
            long long int m_exchange_interval;

            ofstream m_swapfile; // Only used by master

            // Could be safer to have only one vector instead of seperating 
            // control and dependent. Then no chance of using index on wrong vector
            // Vectors of current replica quantities
            vector<double> m_replica_control_qs = vector<double>(3);
            vector<double> m_replica_dependent_qs = vector<double>(3);

            // Vectors of quantities accross all replicas (only filled by master)
            vector<vector<double>> m_control_qs {};

            // Index into the control qs to replica with those qs
            vector<int> m_q_to_repi;

            // Indices into control quantities vector for type
            int m_temp_i {0};
            int m_staple_u_i {1};
            int m_volume_i {2};
            int m_bias_mult_i {3};

            // Indices into dependent quantities vector for type
            int m_enthalpy_i {0};
            int m_staples_i {1};
            int m_bias_i {2};

            // Indices of quantities that will be exchanged
            vector<int> m_exchange_q_is;

            // Initialization methods
            void initialize_control_qs(InputParameters& params);
            void initialize_swap_file(InputParameters& params);

            // Communication methods
            void slave_send(int swap_i);
            void slave_receive(int swap_i);
            void master_receive(int swap_i, vector<vector<double>>& quantities);
            void master_send(int swap_i);
            void master_get_dependent_qs(vector<vector<double>>&);
            virtual void update_control_qs() = 0;
            void update_dependent_qs();

            // Exchange methods
            void attempt_exchange(
                    int swap_i,
                    vector<int>& attempt_count,
                    vector<int>& swap_count);
            bool test_acceptance(double acceptance_p);
            double calc_acceptance_p(
                    vector<pair<double, double>> control_q_pairs,
                    vector<pair<double, double>> dependent_q_pairs); // rep1 and rep2 values

            // Output methods
            void write_swap_entry();
            void write_acceptance_freqs(
                    vector<int> attempt_count,
                    vector<int> swap_count);

            void update_internal(long long int) {};
    };

    class TPTGCMCSimulation: public PTGCMCSimulation {
        // Exchange temperatures
        public:
            TPTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);

        private:
            void update_control_qs();
    };

    class UTPTGCMCSimulation: public PTGCMCSimulation {
        // Exchange temperatures and staple chemical potentials
        public:
            UTPTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);

        private:
            void update_control_qs();
    };

    class HUTPTGCMCSimulation: public PTGCMCSimulation {
        // Exchange temperatures, staple chemical potentials, and bias multpliers
        public:
            HUTPTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);

        private:
            void initialize_exchange_vector();
            void update_control_qs();
    };
}

#endif // PTMC_SIMULATION_H
