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

#include "ceres/ceres.h"

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
            virtual void update_internal() = 0;

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
            int m_exchange_interval;

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
            int m_bias_mult_i {2};

            // Indices into dependent quantities vector for type
            int m_enthalpy_i {0};
            int m_staples_i {1};
            int m_bias_i {2};

            // Indices of quantities that will be exchanged
            vector<int> m_exchange_q_is;

            // Initialization methods
            void initialize_control_qs(InputParameters& params);

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

            void update_internal() {};
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

    using GridPoint = vector<int>;
    using SetOfGridPoints = set<GridPoint>;
    using ArrayOfSets = vector<SetOfGridPoints>;
    using GridInts = unordered_map<GridPoint, int>;
    using GridFloats = unordered_map<GridPoint, double>;
    using GridOfIntArrays = unordered_map<GridPoint, vector<int>>;
    using GridOfFloatArrays = unordered_map<GridPoint, vector<double>>;

    class UmbrellaSamplingSimulation: public GCMCSimulation {
        // For 2D order parameters only
        public:
            UmbrellaSamplingSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            ~UmbrellaSamplingSimulation();
            void run();

        private:
            InputParameters& m_params;
            int m_num_iters;
            long int m_steps;
            SystemOrderParams* m_system_order_params;
            vector<OrderParam*> m_grid_params {};
            GridBiasFunction* m_grid_bias;
            vector<int> m_equil_dif;
            double m_max_D_bias;

            ofstream m_solver_file;
            
            // Variable names set to be consistent with mezei1987
            SetOfGridPoints m_S_n {}; // all grid points visited up to iteration n
            ArrayOfSets m_s_n {}; // grid points visited at iteration n
            SetOfGridPoints m_s_i {}; // grid points visited at current iteration
            GridOfIntArrays m_f_n {}; // number of visits at each grid point at each iteration
            GridInts m_f_i {}; // number of visits at each grid point at current iteration
            GridOfIntArrays m_F_n {}; // total visits at each grid point over all iterations
            GridOfFloatArrays m_r_n {}; // relative contribution of grid points across iterations
            GridOfFloatArrays m_w_n {}; // relative contribution of grid points for each iteration

            GridOfFloatArrays m_p_n {}; // unnormalized weight of grid point for each iteration
            GridFloats m_P_n {}; // unnormalized weight of grid point for all iterations
            vector<double> m_N {}; // normalization of each each iteration
            GridFloats m_E_w {}; // biases at each each iteration

            void update_internal();
            void estimate_current_weights(
                    int n,
                    vector<GridPoint> new_points,
                    vector<GridPoint> old_only_points);
            bool iteration_equilibrium_step(
                    vector<GridPoint> new_points,
                    vector<GridPoint> old_points);
            void estimate_normalizations(int n);
            double estimate_initial_normalization(int n);
            void estimate_final_normalizations();
            void update_grids(
                    int n,
                    vector<GridPoint> new_points,
                    vector<GridPoint> old_points,
                    vector<GridPoint> old_only_points);
            void fill_grid_sets(
                    vector<GridPoint>& new_points,
                    vector<GridPoint>& old_points,
                    vector<GridPoint>& old_only_points);
            void update_bias();
            void output_summary(int n);
    };
    
    class RDevSquareSumCostFunction: public ceres::CostFunction {
        public: 
            RDevSquareSumCostFunction(int num_residuals, vector<int>& parameter_block_sizes);
            ~RDevSquareSumCostFunction() {}
            bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const;
    };

    struct RDevSquareSum {
        // For solving the relative deviation square sum with ceres
        template <typename T> 
        bool operator() (
                T const* const* params,
                T* residual) const;
    };

    GridPoint find_closest_point(
            set<GridPoint> search_set,
            GridPoint target_point,
            int dim);
}

#endif // SIMULATION_H
