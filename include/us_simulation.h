// us_simulation.h

#ifndef US_SIMULATION_H
#define US_SIMULATION_H

#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "boost/mpi/communicator.hpp"
#include "boost/mpi/environment.hpp"
#include "boost/serialization/vector.hpp"

#include "order_params.h"
#include "origami_system.h"
#include "parser.h"
#include "simulation.h"

namespace us {

using std::cout;
using std::ofstream;
using std::ostream;
using std::pair;
using std::set;
using std::string;
using std::unordered_map;
using std::vector;

namespace mpi = boost::mpi;

using biasFunctions::GridBiasFunction;
using biasFunctions::SystemBiases;
using origami::Chain;
using origami::Chains;
using orderParams::OrderParam;
using orderParams::SystemOrderParams;
using origami::OrigamiSystem;
using parser::InputParameters;
using simulation::GCMCSimulation;

using GridPoint = vector<int>;
using SetOfGridPoints = set<GridPoint>;
using ArrayOfSets = vector<SetOfGridPoints>;
using GridInts = unordered_map<GridPoint, int>;
using GridFloats = unordered_map<GridPoint, double>;
using GridOfIntArrays = unordered_map<GridPoint, vector<int>>;
using GridOfFloatArrays = unordered_map<GridPoint, vector<double>>;

// Adaptive US base class. For 2D order parameters only
class USGCMCSimulation: public GCMCSimulation {
  public:
    USGCMCSimulation(
            OrigamiSystem& origami_system,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params);
    virtual ~USGCMCSimulation() {}
    void run();
    void run_equilibration();
    void prepare_iteration(int n);
    void run_simulation(long long int steps);
    void process_iteration(int n);
    bool weights_converged();
    void prepare_production(int n);
    void process_production(int n);
    vector<GridPoint> get_points();
    GridPoint get_current_point();

    // Probably move these to files and interface with grid bias
    void read_weights(string filename);
    void output_weights();

    void set_config_from_traj(string filename, int step);
    void set_config_from_chains(Chains chains);
    void set_output_stream(ostream* out_stream);
    int get_grid_dim();

    long long int m_steps;
    GridBiasFunction& m_grid_bias;
    GridFloats m_E_w {}; // biases at each each iteration

  protected:
    InputParameters& m_params;
    long int m_max_num_iters;
    long long int m_equil_steps;
    long long int m_max_equil_dur;
    long long int m_iter_steps;
    long long int m_max_iter_dur;
    long long int m_prod_steps;
    long long int m_max_prod_dur;
    double m_max_D_bias;
    vector<int> m_equil_dif;

    ostream* m_us_stream {&cout};

    vector<GridPoint> m_new_points {};
    vector<GridPoint> m_old_points {};
    vector<GridPoint> m_old_only_points {};

    SetOfGridPoints m_S_n {}; // all grid points visited up to iteration n
    SetOfGridPoints m_s_i {}; // grid points visited at current iteration
    GridInts m_f_i {}; // number of visits at each grid point at current
                       // iteration
    GridFloats m_p_i {}; // locally normalized weight of grid point for current
                         // iteration
    GridFloats m_old_p_i {}; // previous m_p_i
    GridFloats m_w_i {}; // relative contribution of grid points for current
                         // iteration

    vector<GridPoint> m_points; // timeseries of grid points visited

    void clear_grids();
    void update_internal(long long int step);
    void estimate_current_weights();
    virtual void update_grids(int n) = 0;
    void fill_grid_sets();
    virtual void update_bias(int n) = 0;
    virtual void output_summary(int n) = 0;
};

// A simple adaptive US. Could probably merge into above
class SimpleUSGCMCSimulation: public USGCMCSimulation {
  public:
    SimpleUSGCMCSimulation(
            OrigamiSystem& origami,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params);
    ~SimpleUSGCMCSimulation() {}

  private:
    void update_bias(int);
    void update_grids(int);
    void output_summary(int n);
};

// This doesn't really need to inherit the whole GCSimulation
// Consider a base interface class Simulation
class MWUSGCMCSimulation: public GCMCSimulation {
  public:
    MWUSGCMCSimulation(
            OrigamiSystem& origami_system,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params);
    ~MWUSGCMCSimulation();
    virtual void run();

  protected:
    void setup_window_variables();
    void setup_window_restraints();
    void setup_window_sims(OrigamiSystem& origami);
    void output_iter_summary(int n);
    void update_internal(long long int) {}
    void parse_windows_file(string filename);
    void update_master_order_params(int n);
    void update_master_converged_sims(bool sim_converged, int n);
    void copy_files_to_central_dir(int n);
    void update_starting_config(int n);
    void select_starting_configs(int n);
    void sort_configs_by_ops();

    // MPI variables
    mpi::environment m_env;
    mpi::communicator m_world;
    int m_rank {m_world.rank()};
    int m_master_node {0};

    InputParameters& m_params;
    long int m_max_num_iters;
    string m_local_dir;
    string m_central_dir;

    ofstream* m_us_stream;

    USGCMCSimulation* m_us_sim;
    int m_windows {0};
    int m_grid_dim {0};
    vector<string> m_window_bias_tags {};
    vector<GridPoint> m_window_mins {};
    vector<GridPoint> m_window_maxs {};
    vector<unsigned int> m_num_points {}; // num grid points per window
    vector<string> m_output_filebases {};
    vector<string> m_window_postfixes {};
    string m_starting_file {};
    int m_starting_step {};

    vector<vector<GridPoint>> m_points {};
    unordered_map<GridPoint, vector<pair<int, int>>>
            m_order_param_to_configs {};
    vector<bool> m_sims_converged {};
    int m_num_sims_converged {0};
    vector<int> m_current_iters {};
    vector<string> m_starting_files {};
    vector<int> m_starting_steps {};
};

class PTMWUSGCMCSimulation: public MWUSGCMCSimulation {
  public:
    PTMWUSGCMCSimulation(
            OrigamiSystem& origami_system,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params);
    void run() override;

  private:
    void run_swaps(long long int swaps, long long int dur);
    void slave_send_ops(int swap_i);
    bool slave_send_and_recieve_chains(int swap_i);
    void attempt_exchange(int swap_i);
    void master_send_kill(int swap_i);
    void write_acceptance_freqs();
    void initialize_swap_file(InputParameters& params);
    void write_swap_entry(long long int step);
    bool test_acceptance(double p_accept);

    long long int m_iter_swaps;
    long long int m_max_iter_dur;
    long long int m_prod_swaps;
    long long int m_max_prod_dur;
    long long int m_exchange_interval;
    vector<int> m_attempt_count;
    vector<int> m_swap_count;
    vector<int> m_win_to_configi;
    ofstream m_swapfile; // Only used by master
};
} // namespace us

#endif // US_SIMULATION_H
