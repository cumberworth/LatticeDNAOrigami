// simulation.h

#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/process.hpp>

#include "bias_functions.h"
#include "files.h"
#include "ideal_random_walk.h"
#include "movetypes.h"
#include "order_params.h"
#include "origami_system.h"
#include "parser.h"
#include "random_gens.h"

namespace simulation {

using std::map;
using std::ofstream;
using std::ostream;
using std::set;
using std::string;
using std::unique_ptr;
using std::unordered_map;
using std::vector;
using std::chrono::steady_clock;

namespace bp = boost::process;

using biasFunctions::SystemBiases;
using origami::Chain;
using origami::Chains;
using files::OrigamiMovetypeFile;
using files::OrigamiOrientationOutputFile;
using files::OrigamiOutputFile;
using files::OrigamiStateOutputFile;
using files::OrigamiVCFOutputFile;
using files::OrigamiVSFOutputFile;
using idealRandomWalk::IdealRandomWalks;
using movetypes::MCMovetype;
using orderParams::SystemOrderParams;
using origami::OrigamiSystem;
using parser::InputParameters;
using randomGen::RandomGens;

vector<OrigamiOutputFile*> setup_output_files(
        InputParameters& params,
        string output_filebase,
        OrigamiSystem& origami,
        SystemOrderParams& ops,
        SystemBiases& biases,
        RandomGens& random_gens);

void setup_config_files(
        const string filebase,
        const int max_total_staples,
        const int max_staple_size,
        const int freq,
        OrigamiSystem& origami,
        vector<OrigamiOutputFile*>& files);

class GCMCSimulation {
  public:
    GCMCSimulation(
            OrigamiSystem& origami_system,
            SystemOrderParams& ops,
            SystemBiases& biases,
            InputParameters& params);
    virtual ~GCMCSimulation();
    virtual void run() = 0;
    Chains get_chains();
    long long int simulate(
            long long int steps,
            long long int start_step = 0,
            bool summarize = true,
            steady_clock::time_point = steady_clock::now());

  protected:
    // Shared interface
    virtual void update_internal(long long int step) = 0;

    // Shared methods
    void construct_movetypes(InputParameters& params);
    MCMovetype* setup_orientation_movetype(
            int i,
            string type,
            string label,
            OrigamiMovetypeFile& movetypes_file,
            MCMovetype* movetype);
    MCMovetype* setup_staple_exchange_movetype(
            int i,
            string type,
            string label,
            OrigamiMovetypeFile& movetypes_file,
            MCMovetype* movetype);
    MCMovetype* setup_staple_regrowth_movetype(
            int i,
            string type,
            string label,
            OrigamiMovetypeFile& movetypes_file,
            MCMovetype* movetype);
    MCMovetype* setup_scaffold_regrowth_movetype(
            int i,
            string type,
            string label,
            OrigamiMovetypeFile& movetypes_file,
            MCMovetype* movetype);
    MCMovetype* setup_scaffold_transform_movetype(
            int i,
            string type,
            string label,
            OrigamiMovetypeFile& movetypes_file,
            MCMovetype* movetype);
    MCMovetype* setup_rg_movetype(
            int i,
            string type,
            string label,
            OrigamiMovetypeFile& movetypes_file,
            MCMovetype* movetype);
    void set_max_dur(long long int dur);
    MCMovetype& select_movetype();
    void write_log_entry(
            const long long int step,
            bool accepted,
            MCMovetype& movetype);
    void write_log_summary();
    void setup_vmd_pipe();
    void pipe_to_vmd();
    void close_vmd_pipe();
    void close_output_files();

    OrigamiSystem& m_origami_system;
    SystemOrderParams& m_ops;
    SystemBiases& m_biases;

    ostream* m_logging_stream;
    int m_logging_freq;
    int m_centering_freq;
    int m_centering_domain;
    int m_constraint_check_freq;
    int m_vmd_pipe_freq;
    double m_max_duration;
    InputParameters& m_params;
    vector<OrigamiOutputFile*> m_output_files {};
    vector<OrigamiOutputFile*> m_config_per_move_files {};
    vector<unique_ptr<MCMovetype>> m_movetypes {};
    vector<double> m_movetype_freqs {};
    vector<double> m_cumulative_probs {};
    RandomGens m_random_gens {};
    IdealRandomWalks m_ideal_random_walks {};

    // VMD realtime visualization
    OrigamiVSFOutputFile* vmd_struct_file {nullptr};
    OrigamiVCFOutputFile* vmd_coors_file {nullptr};
    OrigamiStateOutputFile* vmd_states_file {nullptr};
    OrigamiOrientationOutputFile* vmd_ores_file {nullptr};
    bp::child* vmd_proc {nullptr};
};
} // namespace simulation

#endif // SIMULATION_H
