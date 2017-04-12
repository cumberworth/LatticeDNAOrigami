// us_simulation.h

#ifndef US_SIMULATION_H
#define US_SIMULATION_H

#include "boost/serialization/vector.hpp"

#include "simulation.h"

using namespace Simulation;

using GridPoint = vector<int>;
using SetOfGridPoints = set<GridPoint>;
using ArrayOfSets = vector<SetOfGridPoints>;
using GridInts = unordered_map<GridPoint, int>;
using GridFloats = unordered_map<GridPoint, double>;
using GridOfIntArrays = unordered_map<GridPoint, vector<int>>;
using GridOfFloatArrays = unordered_map<GridPoint, vector<double>>;

namespace US {

    GridPoint find_closest_point(
            set<GridPoint> search_set,
            GridPoint target_point,
            int dim);

    class USGCMCSimulation: public GCMCSimulation {
        // Adaptive US base class. For 2D order parameters only
        public:
            USGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            virtual ~USGCMCSimulation();
            void run();
            void run_equilibration();
            void run_iteration(int n);
            bool weights_converged();
            void run_production(int n);
            vector<GridPoint> get_points();
            void read_weights(string filename);
            void output_weights();
            void set_config_from_traj(string filename, int step);
            void set_output_stream(ostream* out_stream);
            int get_grid_dim();

        protected:
            InputParameters& m_params;
            long int m_max_num_iters;
            long int m_equil_steps;
            long int m_steps;
            long int m_prod_steps;
            double m_max_rel_P_diff;
            SystemOrderParams* m_system_order_params;
            vector<OrderParam*> m_grid_params {};
            GridBiasFunction* m_grid_bias;
            double m_max_D_bias;
            vector<int> m_equil_dif;

            ostream* m_us_stream {&cout};

            vector<GridPoint> m_new_points {};
            vector<GridPoint> m_old_points {};
            vector<GridPoint> m_old_only_points {};

            // Variable names set to be consistent with mezei1987
            SetOfGridPoints m_S_n {}; // all grid points visited up to iteration n
            SetOfGridPoints m_s_i {}; // grid points visited at current iteration
            GridInts m_f_i {}; // number of visits at each grid point at current iteration
            GridFloats m_p_i {}; // locally normalized weight of grid point for current iteration
            GridFloats m_lP_n {}; // locally normalized weight of grid point
            GridFloats m_old_lP_n {}; // previous m_lP_n
            GridFloats m_w_i {}; // relative contribution of grid points for current iteration
            GridFloats m_E_w {}; // biases at each each iteration

            vector<GridPoint> m_points; // timeseries of grid points visited

            void clear_grids();
            void update_internal(long int step);
            void estimate_current_weights();
            virtual void update_grids(int n) = 0;
            bool iteration_equilibrium_step();
            void fill_grid_sets();
            virtual void update_bias(int n) = 0;
            virtual void output_summary(int n) = 0;
            void prod_output_summary();
    };

    class SimpleUSGCMCSimulation: public USGCMCSimulation {
        // Simple adaptive US
        public:
            SimpleUSGCMCSimulation(
                    OrigamiSystem& origami,
                    InputParameters& params);
            ~SimpleUSGCMCSimulation() {}
            
        private:
            void update_bias(int);
            void update_grids(int);
            void output_summary(int n);
    };

    class MWUSGCMCSimulation: public GCMCSimulation {
        public:
            MWUSGCMCSimulation(
                    OrigamiSystem& origami_system,
                    InputParameters& params);
            ~MWUSGCMCSimulation();
            void run();

        private:

            // MPI variables
            mpi::environment m_env;
            mpi::communicator m_world;
            int m_rank {m_world.rank()};
            int m_master_node {0};

            InputParameters& m_params;
            long int m_max_num_iters;
            long int m_equil_steps;
            long int m_steps;
            long int m_prod_steps;
            SystemBiases* m_system_biases;

            ofstream* m_us_stream;

            USGCMCSimulation* m_us_sim;
            int m_windows {0};
            int m_grid_dim {0};
            vector<GridPoint> m_window_mins {};
            vector<GridPoint> m_window_maxs {};
            vector<string> m_output_filebases {};
            vector<string> m_window_postfixes {};
            string m_starting_file {};
            int m_starting_step {};

            vector<vector<GridPoint>> m_points {};
            unordered_map<GridPoint, vector<pair<int, int>>> m_order_param_to_configs {};
            vector<bool> m_sims_converged {};
            int m_num_sims_converged {0};
            vector<int> m_current_iters {};
            vector<string> m_starting_files {};
            vector<int> m_starting_steps {};

            void update_internal(long int) {}
            void parse_windows_file(string filename);
            void update_master_order_params(int n);
            void update_master_converged_sims(bool sim_converged, int n);
            void update_starting_config(int n);
            void select_starting_configs();
    };
}

#endif // US_SIMULATION_H
