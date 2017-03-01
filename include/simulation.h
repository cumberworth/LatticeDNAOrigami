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

#include "parser.h"
#include "origami_system.h"
#include "movetypes.h"
#include "files.h"
#include "ideal_random_walk.h"
#include "order_params.h"

using std::vector;
using std::unique_ptr;
using std::map;
using std::ostream;
using std::unordered_map;
using std::set;

namespace mpi = boost::mpi;

using namespace Parser;
using namespace Origami;
using namespace Movetypes;
using namespace Files;
using namespace IdealRandomWalk;
using namespace OrderParams;

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

            // Shared interface
            virtual void update_internal();

            // Shared methods
            void simulate(long int steps, int start_step=0);
            unique_ptr<MCMovetype> select_movetype();
            void write_log_entry(
                    long int step,
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
            void update_internal() {};
            long int m_steps;
    };

    class AnnealingGCMCSimulation: public GCMCSimulation {
        public:
            AnnealingGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            void run();
        private:
            void update_internal() {};
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
            ~PTGCMCSimulation();
            void run();

        protected:
            mpi::environment m_env;
            mpi::communicator m_world;
            int m_rank {m_world.rank()};
            int m_master_rep {0};
            int m_swaps;
            int m_num_reps;
            int m_exchange_interval;
            vector<vector<double>> m_control_qs;
            vector<vector<double>> m_dependent_qs;

            void update_internal() {};

            // Data only filled in master
            ofstream m_swapfile;
            vector<int> m_q_to_repi; // Index in quantity to swap to replica number with that value

            virtual void slave_send(int swap_i) = 0;
            virtual void slave_recieve(int swap_i) = 0;
            void master_receive(int swap_i, vector<vector<double>>& quantities);
            void master_send(int swap_i);
            virtual void master_set_control_qs() = 0;
            virtual void master_get_dependent_qs(vector<vector<double>>&) = 0;
            void attempt_exchange(
                    int swap_i,
                    vector<int>& attempt_count,
                    vector<int>& swap_count);
            virtual void update_control_qs() = 0;
            bool test_acceptance(double acceptance_p);
            virtual double calc_acceptance_p(
                    vector<pair<double, double>> control_q_pairs,
                    vector<pair<double, double>> dependent_q_pairs) = 0; // rep1 and rep2 values
            void write_swap_entry();
            void write_acceptance_freqs(
                    vector<int> attempt_count,
                    vector<int> swap_count);
    };

    class TPTGCMCSimulation: public PTGCMCSimulation {
        public:
            TPTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);

        private:
            double m_temp;

            // Indices into quantities vector for type
            int m_temp_i {0};

            // Indices into dependent quantities vector for type
            int m_enthalpy_i {0};

            void slave_send(int swap_i);
            void master_set_control_qs();
            void update_control_qs();
            double calc_acceptance_p(
                    vector<pair<double, double>> control_q_pairs,
                    vector<pair<double, double>> dependent_q_pairs);
    };

    class UTPTGCMCSimulation: public PTGCMCSimulation {
        public:
            UTPTGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);

        private:
            double m_temp;
            double m_staple_u;

            // Indices into control quantities vector for type
            int m_temp_i {0};
            int m_staple_u_i {1};

            // Indices into dependent quantities vector for type
            int m_enthalpy_i {0};
            int m_staples_i {1};

            void slave_send(int swap_i);
            void master_set_control_qs();
            void update_control_qs();
            double calc_acceptance_p(
                    vector<pair<double, double>> control_q_pairs,
                    vector<pair<double, double>> dependent_q_pairs);
    };
            vector<double> m_bias_mults;
            double m_volume;
            double m_bias_mult;
            double m_staple_u;

    using GridPoint = vector<int>;
    using SetOfGridPoints = set<GridPoint>;
    using GridInts = unordered_map<GridPoint, int>;
    using GridFloats = unordered_map<GridPoint, double>;
    using GridOfIntArrays = unordered_map<GridPoint, vector<int>>;
    using GridOfFloatArrays = unordered_map<GridPoint, vector<double>>;
    using ArrayOfSets = vector<SetOfGridPoints>;

    class UmbrellaSamplingSimulation: public GCMCSimulation {
        // For 2D order parameters only at this point
        public:
            UmbrellaSamplingSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            ~UmbrellaSamplingSimulation();
            void run();

        private:
            SystemOrderParams m_system_order_params;
            SystemBiases m_system_biases;

            ofstream m_solver_file;
            
            // Variable names set to be consistent with mezei1987
            SetOfGridPoints m_S_n {}; // grid points visited in complete interations
            ArrayOfSets m_s_n {}; // grid points visisted at iteration n
            GridOfIntArrays m_f_n {}; // number of visits at each grid point at each iteration
            GridInts m_F_n {}; // total visits at each grid point over all iterations
            vector<double> m_N {}; // normalization of each grid point, for each iteration
            GridOfFloatArrays m_r_n {}; // relative contribution of grid points across itertions
            GridOfFloatArrays m_w_n {}; // relative contribution of grid points for each iteration
            GridOfFloatArrays m_p_n {}; // unnormalized weight of grid point for each iteration
            GridInts m_P_n {}; // unnormalized weight of grid point for all iterations
            GridOfFloatArrays m_E_w {}; // biases at each each iteration

            void update_internal() {};
            void estimate_current_weights();
            bool iteration_equilibrium_step();
            void estimate_normalizations(int n);
            double estimate_initial_normalization(int n);
            void estimate_final_normalizations();
    };

    struct RDevSquareSum {
        template <typename T> 
        bool operator() (
                T const* const* parameters,
                T* residual) const;
    };

    int find_closest_point(set<GridPoint> search_set, int dim);
}

#endif // SIMULATION_H
