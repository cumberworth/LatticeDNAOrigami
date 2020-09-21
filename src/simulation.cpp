// simulation.cpp

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <utility>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/process.hpp>

#include "cb_movetypes.h"
#include "files.h"
#include "met_movetypes.h"
#include "movetypes.h"
#include "orientation_movetype.h"
#include "random_gens.h"
#include "rg_movetypes.h"
#include "simulation.h"
#include "transform_movetypes.h"
#include "utility.h"

#include "random_gens.h"

namespace simulation {

using std::cout;
using std::min;
using std::pair;
using std::pow;
using std::setw;
using std::chrono::steady_clock;

namespace mpi = boost::mpi;
namespace bp = boost::process;

using files::OrigamiCountsOutputFile;
using files::OrigamiEnergiesOutputFile;
using files::OrigamiMovetypeFile;
using files::OrigamiOrderParamsOutputFile;
using files::OrigamiStaplesBoundOutputFile;
using files::OrigamiStaplesFullyBoundOutputFile;
using files::OrigamiTimesOutputFile;
using files::OrigamiTrajOutputFile;
using files::RandomEngineStateInputFile;
using files::RandomEngineStateOutputFile;

vector<OrigamiOutputFile*> setup_output_files(
        InputParameters& params,
        string output_filebase,
        OrigamiSystem& origami,
        SystemOrderParams& ops,
        SystemBiases& biases,
        RandomGens& random_gens) {

    // Hack to get a vsf file
    OrigamiVSFOutputFile vsf_file {
            output_filebase + ".vsf",
            0,
            params.m_max_total_staples,
            params.m_max_staple_size,
            origami};
    vsf_file.write(0, 0);

    vector<OrigamiOutputFile*> outs {};
    if (params.m_configs_output_freq != 0) {
        OrigamiOutputFile* config_out = new OrigamiTrajOutputFile {
                output_filebase + ".trj",
                params.m_configs_output_freq,
                params.m_max_total_staples,
                params.m_max_staple_size,
                origami};
        outs.push_back(config_out);
    }
    if (params.m_vtf_output_freq != 0) {
        setup_config_files(
                output_filebase,
                params.m_max_total_staples,
                params.m_max_staple_size,
                params.m_vtf_output_freq,
                origami,
                outs);
    }
    if (params.m_counts_output_freq != 0) {
        OrigamiOutputFile* counts_out = new OrigamiCountsOutputFile {
                output_filebase + ".counts",
                params.m_counts_output_freq,
                params.m_max_total_staples,
                params.m_max_staple_size,
                origami};
        outs.push_back(counts_out);
        OrigamiOutputFile* staples_bound = new OrigamiStaplesBoundOutputFile {
                output_filebase + ".staples",
                params.m_counts_output_freq,
                params.m_max_total_staples,
                params.m_max_staple_size,
                origami};
        outs.push_back(staples_bound);
        OrigamiOutputFile* staples_fully_bound =
                new OrigamiStaplesFullyBoundOutputFile {
                        output_filebase + ".staplestates",
                        params.m_counts_output_freq,
                        params.m_max_total_staples,
                        params.m_max_staple_size,
                        origami};
        outs.push_back(staples_fully_bound);
    }
    if (params.m_times_output_freq != 0) {
        OrigamiOutputFile* times_out = new OrigamiTimesOutputFile {
                output_filebase + ".times",
                params.m_times_output_freq,
                params.m_max_total_staples,
                params.m_max_staple_size,
                origami};
        outs.push_back(times_out);
    }
    if (params.m_energies_output_freq != 0) {
        OrigamiOutputFile* energies_out = new OrigamiEnergiesOutputFile {
                output_filebase + ".ene",
                params.m_energies_output_freq,
                params.m_max_total_staples,
                params.m_max_staple_size,
                origami,
                biases};
        outs.push_back(energies_out);
    }
    if (params.m_order_params_output_freq != 0) {
        OrigamiOutputFile* order_params_out = new OrigamiOrderParamsOutputFile {
                output_filebase + ".ops",
                params.m_order_params_output_freq,
                params.m_max_total_staples,
                params.m_max_staple_size,
                origami,
                ops,
                params.m_ops_to_output};
        outs.push_back(order_params_out);
    }
    if (params.m_rand_engine_state_output_freq != 0) {
        OrigamiOutputFile* rand_state_out = new RandomEngineStateOutputFile {
                output_filebase + ".randstate",
                params.m_rand_engine_state_output_freq,
                random_gens,
                origami};
        outs.push_back(rand_state_out);
    }

    return outs;
}

void setup_config_files(
        const string filebase,
        const int max_total_staples,
        const int max_staple_size,
        const int freq,
        OrigamiSystem& origami,
        vector<OrigamiOutputFile*>& files) {

    OrigamiOutputFile* config_out = new OrigamiVCFOutputFile {
            filebase + ".vcf",
            freq,
            max_total_staples,
            max_staple_size,
            origami};
    files.push_back(config_out);

    config_out = new OrigamiStateOutputFile {
            filebase + ".states",
            freq,
            max_total_staples,
            max_staple_size,
            origami};
    files.push_back(config_out);

    config_out = new OrigamiOrientationOutputFile {
            filebase + ".ores",
            freq,
            max_total_staples,
            max_staple_size,
            origami};
    files.push_back(config_out);
}

GCMCSimulation::GCMCSimulation(
        OrigamiSystem& origami_system,
        SystemOrderParams& ops,
        SystemBiases& biases,
        InputParameters& params):
        m_origami_system {origami_system},
        m_ops {ops},
        m_biases {biases},
        m_params {params} {

    m_logging_freq = params.m_logging_freq;
    m_centering_freq = params.m_centering_freq;
    m_centering_domain = params.m_centering_domain;
    m_constraint_check_freq = params.m_constraint_check_freq;
    m_vmd_pipe_freq = params.m_vmd_pipe_freq;
    m_max_duration = params.m_max_duration;

    // This only makes sense for serial debugging
    if (m_params.m_random_seed != -1) {
        m_random_gens.set_seed(m_params.m_random_seed);
        cout << "Using specified seed: " << m_params.m_random_seed << "\n";
    }
    else if (not m_params.m_rand_engine_state_file.empty()) {
        cout << "Loading random engine state\n";
        RandomEngineStateInputFile rand_engine_state_file {
                m_params.m_rand_engine_state_file};
        std::istringstream rand_engine_state {
                rand_engine_state_file.read_state(m_params.m_restart_step)};
        rand_engine_state.setf(std::ios::dec);
        rand_engine_state >> m_random_gens.m_random_engine;
    }

    if (m_vmd_pipe_freq != 0) {
        setup_vmd_pipe();
    }

    // HACK (files won't be right if filebase changed in derived constructor)
    if (params.m_vcf_per_domain) {
        setup_config_files(
                params.m_output_filebase + "_move",
                params.m_max_total_staples,
                params.m_max_staple_size,
                params.m_configs_output_freq,
                m_origami_system,
                m_config_per_move_files);
    }

    // Constructor movetypes
    construct_movetypes(params);

    // Create cumulative probability array
    double cum_prob {0};
    for (size_t i {0}; i != m_movetypes.size(); i++) {
        cum_prob += m_movetype_freqs[i];
        m_cumulative_probs.push_back(cum_prob);
    }

    // Load precalculated ideal random walk count data
    if (params.m_num_walks_filename.size() != 0) {
        std::ifstream num_walks_file {params.m_num_walks_filename};
        boost::archive::binary_iarchive num_walks_arch {num_walks_file};
        num_walks_arch >> m_ideal_random_walks;
    }
}

GCMCSimulation::~GCMCSimulation() {
    close_output_files();
    if (m_vmd_pipe_freq != 0) {
        close_vmd_pipe();
    }
}

void GCMCSimulation::construct_movetypes(InputParameters& params) {
    OrigamiMovetypeFile movetypes_file {params.m_movetype_filename};
    vector<string> types {movetypes_file.get_types()};
    vector<string> labels {movetypes_file.get_labels()};
    m_movetype_freqs = movetypes_file.get_freqs();
    for (size_t i {0}; i != types.size(); i++) {
        MCMovetype* movetype {nullptr};
        string label {labels[i]};

        // This is still sort of ugly. What if I end up wanting to share
        // options between all CB moves, or between scaffold regrowth and
        // scaffold transform, for example?
        string type {types[i]};
        if (type == "OrientationRotation") {
            movetype = setup_orientation_movetype(
                    i, type, label, movetypes_file, movetype);
        }
        else if (type == "MetStapleExchange") {
            movetype = setup_staple_exchange_movetype(
                    i, type, label, movetypes_file, movetype);
        }
        else if (type == "MetStapleRegrowth" or type == "CBStapleRegrowth") {
            movetype = setup_staple_regrowth_movetype(
                    i, type, label, movetypes_file, movetype);
        }
        else if (
                type == "CTCBScaffoldRegrowth" or
                type == "CTCBJumpScaffoldRegrowth") {
            movetype = setup_scaffold_regrowth_movetype(
                    i, type, label, movetypes_file, movetype);
        }
        else if (
                type == "CTCBLinkerRegrowth" or
                type == "CTCBClusteredLinkerRegrowth" or
                type == "CTRGLinkerRegrowth" or
                type == "CTRGClusteredLinkerRegrowth") {

            movetype = setup_scaffold_transform_movetype(
                    i, type, label, movetypes_file, movetype);
        }
        else if (
                type == "CTRGScaffoldRegrowth" or
                type == "CTRGJumpScaffoldRegrowth") {
            movetype =
                    setup_rg_movetype(i, type, label, movetypes_file, movetype);
        }
        else {
            cout << "No such movetype\n";
            throw utility::SimulationMisuse {};
        }
        m_movetypes.emplace_back(movetype);
    }
}

void GCMCSimulation::set_max_dur(long long int dur) { m_max_duration = dur; }

MCMovetype* GCMCSimulation::setup_orientation_movetype(
        int,
        string,
        string label,
        OrigamiMovetypeFile&,
        MCMovetype* movetype) {

    movetype = new movetypes::OrientationRotationMCMovetype {
            m_origami_system,
            m_random_gens,
            m_ideal_random_walks,
            m_config_per_move_files,
            label,
            m_ops,
            m_biases,
            m_params};

    return movetype;
}

MCMovetype* GCMCSimulation::setup_staple_exchange_movetype(
        int i,
        string,
        string label,
        OrigamiMovetypeFile& movetypes_file,
        MCMovetype* movetype) {

    vector<double> exchange_mults {
            movetypes_file.get_double_vector_option(i, "exchange_mults")};
    bool adaptive_exchange {
            movetypes_file.get_bool_option(i, "adaptive_exchange")};
    movetype = new movetypes::MetStapleExchangeMCMovetype {
            m_origami_system,
            m_random_gens,
            m_ideal_random_walks,
            m_config_per_move_files,
            label,
            m_ops,
            m_biases,
            m_params,
            exchange_mults,
            adaptive_exchange};

    return movetype;
}

MCMovetype* GCMCSimulation::setup_staple_regrowth_movetype(
        int,
        string type,
        string label,
        OrigamiMovetypeFile&,
        MCMovetype* movetype) {

    if (type == "MetStapleRegrowth") {
        movetype = new movetypes::MetStapleRegrowthMCMovetype {
                m_origami_system,
                m_random_gens,
                m_ideal_random_walks,
                m_config_per_move_files,
                label,
                m_ops,
                m_biases,
                m_params};
    }
    else if (type == "CBStapleRegrowth") {
        movetype = new movetypes::CBStapleRegrowthMCMovetype {
                m_origami_system,
                m_random_gens,
                m_ideal_random_walks,
                m_config_per_move_files,
                label,
                m_ops,
                m_biases,
                m_params};
    }

    return movetype;
}

MCMovetype* GCMCSimulation::setup_scaffold_regrowth_movetype(
        int i,
        string type,
        string label,
        OrigamiMovetypeFile& movetypes_file,
        MCMovetype* movetype) {

    // int excluded_staples {movetypes_file.get_int_option(i,
    //        "num_excluded_staples")};
    int excluded_staples {0};
    int max_regrowth {movetypes_file.get_int_option(i, "max_regrowth")};
    if (type == "CTCBScaffoldRegrowth") {
        movetype = new movetypes::CTCBScaffoldRegrowthMCMovetype {
                m_origami_system,
                m_random_gens,
                m_ideal_random_walks,
                m_config_per_move_files,
                label,
                m_ops,
                m_biases,
                m_params,
                excluded_staples,
                max_regrowth};
    }
    else if (type == "CTCBJumpScaffoldRegrowth") {
        int max_seg_regrowth {
                movetypes_file.get_int_option(i, "max_seg_regrowth")};
        movetype = new movetypes::CTCBJumpScaffoldRegrowthMCMovetype {
                m_origami_system,
                m_random_gens,
                m_ideal_random_walks,
                m_config_per_move_files,
                label,
                m_ops,
                m_biases,
                m_params,
                excluded_staples,
                max_regrowth,
                max_seg_regrowth};
    }

    return movetype;
}

MCMovetype* GCMCSimulation::setup_scaffold_transform_movetype(
        int i,
        string type,
        string label,
        OrigamiMovetypeFile& movetypes_file,
        MCMovetype* movetype) {

    int excluded_staples {0};
    int max_disp {movetypes_file.get_int_option(i, "max_disp")};
    int max_turns {movetypes_file.get_int_option(i, "max_turns")};
    int max_regrowth {movetypes_file.get_int_option(i, "max_regrowth")};
    int max_linker_length {
            movetypes_file.get_int_option(i, "max_linker_length")};
    int num_transforms {movetypes_file.get_int_option(i, "num_transforms")};
    if (type == "CTCBLinkerRegrowth") {
        movetype = new movetypes::CTCBLinkerRegrowthMCMovetype {
                m_origami_system,
                m_random_gens,
                m_ideal_random_walks,
                m_config_per_move_files,
                label,
                m_ops,
                m_biases,
                m_params,
                excluded_staples,
                max_disp,
                max_turns,
                static_cast<unsigned int>(max_regrowth),
                static_cast<unsigned int>(max_linker_length),
                num_transforms};
    }
    else if (type == "CTCBClusteredLinkerRegrowth") {
        movetype = new movetypes::CTCBClusteredLinkerRegrowthMCMovetype {
                m_origami_system,
                m_random_gens,
                m_ideal_random_walks,
                m_config_per_move_files,
                label,
                m_ops,
                m_biases,
                m_params,
                excluded_staples,
                max_disp,
                max_turns,
                static_cast<unsigned int>(max_linker_length),
                num_transforms};
    }
    else if (type == "CTRGLinkerRegrowth") {
        int max_num_recoils {
                movetypes_file.get_int_option(i, "max_num_recoils")};
        int max_c_attempts {movetypes_file.get_int_option(i, "max_c_attempts")};
        movetype = new movetypes::CTRGLinkerRegrowthMCMovetype {
                m_origami_system,
                m_random_gens,
                m_ideal_random_walks,
                m_config_per_move_files,
                label,
                m_ops,
                m_biases,
                m_params,
                excluded_staples,
                max_num_recoils,
                max_c_attempts,
                max_disp,
                max_turns,
                static_cast<unsigned int>(max_regrowth),
                static_cast<unsigned int>(max_linker_length),
                num_transforms};
    }

    return movetype;
}

MCMovetype* GCMCSimulation::setup_rg_movetype(
        int i,
        string type,
        string label,
        OrigamiMovetypeFile& movetypes_file,
        MCMovetype* movetype) {

    //        int excluded_staples {movetypes_file.get_int_option(i,
    //                "num_excluded_staples")};
    int excluded_staples {0};
    int max_num_recoils {movetypes_file.get_int_option(i, "max_num_recoils")};
    int max_c_attempts {movetypes_file.get_int_option(i, "max_c_attempts")};
    int max_regrowth {movetypes_file.get_int_option(i, "max_regrowth")};
    if (type == "CTRGScaffoldRegrowth") {
        movetype = new movetypes::CTRGScaffoldRegrowthMCMovetype {
                m_origami_system,
                m_random_gens,
                m_ideal_random_walks,
                m_config_per_move_files,
                label,
                m_ops,
                m_biases,
                m_params,
                excluded_staples,
                max_num_recoils,
                max_c_attempts,
                max_regrowth};
    }
    else if (type == "CTRGJumpScaffoldRegrowth") {
        int max_seg_regrowth {
                movetypes_file.get_int_option(i, "max_seg_regrowth")};
        movetype = new movetypes::CTRGJumpScaffoldRegrowthMCMovetype {
                m_origami_system,
                m_random_gens,
                m_ideal_random_walks,
                m_config_per_move_files,
                label,
                m_ops,
                m_biases,
                m_params,
                excluded_staples,
                max_num_recoils,
                max_c_attempts,
                max_regrowth,
                max_seg_regrowth};
    }

    return movetype;
}

long long int GCMCSimulation::simulate(
        long long int steps,
        long long int start_step,
        bool summarize,
        steady_clock::time_point start) {

    long long int step {start_step + 1};
    for (; step != (steps + start_step + 1); step++) {

        // Pick movetype and apply
        MCMovetype& movetype {select_movetype()};
        bool accepted;
        double old_ene {m_origami_system.energy()};
        accepted = movetype.attempt_move(step);
        if (not accepted) {
            // cout << "\nReseting configuration\n";
            movetype.reset_origami();
            m_ops.update_move_params();
            m_biases.calc_move();
        }
        double new_ene {m_origami_system.energy()};
        if (new_ene - old_ene >= 100) {
            cout << "Hang on\n";
        }

        // Center and check constraints
        if (m_centering_freq != 0 and step % m_centering_freq == 0) {
            m_origami_system.center(m_centering_domain);
        }
        if (m_constraint_check_freq != 0 and
            step % m_constraint_check_freq == 0) {
            // cout << "check\n";
            try {
                m_origami_system.check_all_constraints();
            } catch (utility::OrigamiMisuse) {
                write_log_entry(step, accepted, movetype);
                std::chrono::duration<double> dt {
                        (steady_clock::now() - start)};
                for (auto output_file: m_output_files) {
                    output_file->write(step, dt.count());
                }
                pipe_to_vmd();
                cout << "Origami misuse at constraint check\n";
                break;
            }
        }

        std::chrono::duration<double> dt {(steady_clock::now() - start)};
        if (dt.count() > m_max_duration) {
            cout << "Maximum time allowed reached\n";
            break;
        }

        // Write log entry to standard out
        if (m_logging_freq != 0 and step % m_logging_freq == 0) {
            write_log_entry(step, accepted, movetype);
        }

        // VMD pipe
        if (m_vmd_pipe_freq != 0 and step % m_vmd_pipe_freq == 0) {
            pipe_to_vmd();
        }

        // Update internal simulation variables
        update_internal(step);

        // Write system properties to file
        for (auto output_file: m_output_files) {
            if (output_file->m_write_freq != 0 and
                step % output_file->m_write_freq == 0) {
                output_file->write(step, dt.count());
            }
        }
    }
    if (summarize) {
        write_log_summary();
    }

    return step;
}

MCMovetype& GCMCSimulation::select_movetype() {
    double prob {m_random_gens.uniform_real()};
    size_t i;
    for (i = 0; i != m_cumulative_probs.size(); i++) {
        if (prob < m_cumulative_probs[i]) {
            break;
        }
    }

    return *m_movetypes[i];
}

void GCMCSimulation::write_log_entry(
        const long long int step,
        bool accepted,
        MCMovetype& movetype) {

    *m_logging_stream << "Step: " << step << "\n";
    *m_logging_stream << "Temperature: " << m_origami_system.m_temp << "\n";
    *m_logging_stream << "Bound staples: " << m_origami_system.num_staples()
                      << "\n";
    *m_logging_stream << "Unique bound staples: "
                      << m_origami_system.num_unique_staples() << "\n";
    *m_logging_stream << "Fully bound domain pairs: "
                      << m_origami_system.num_fully_bound_domain_pairs()
                      << "\n";
    *m_logging_stream << "Misbound domain pairs: "
                      << m_origami_system.num_misbound_domain_pairs() << "\n";
    *m_logging_stream << "Stacked domain pairs: "
                      << m_origami_system.num_stacked_domain_pairs() << "\n";
    *m_logging_stream << "Linear helix triplets: "
                      << m_origami_system.num_linear_helix_trips() << "\n";
    *m_logging_stream << "Stacked junction quadruplets: "
                      << m_origami_system.num_stacked_junct_quads() << "\n";
    *m_logging_stream << "Staple counts: ";
    for (auto staple_count: m_origami_system.get_staple_counts()) {
        *m_logging_stream << staple_count << " ";
    }
    *m_logging_stream << "\n";
    *m_logging_stream << "System energy: " << m_origami_system.energy() << "\n";
    *m_logging_stream << "Total external bias: " << m_biases.get_total_bias()
                      << "\n";
    *m_logging_stream << "Domain update external bias: "
                      << m_biases.get_domain_update_bias() << "\n";
    *m_logging_stream << "Move update external bias: "
                      << m_biases.get_move_update_bias() << "\n";
    *m_logging_stream << "Movetype: " << movetype.get_label() << "\n";
    *m_logging_stream << "Accepted: " << std::boolalpha << accepted << "\n";
    *m_logging_stream << "\n";
}

void GCMCSimulation::write_log_summary() {
    *m_logging_stream << "Run summary"
                      << "\n\n";
    ofstream* movetype_sum_stream;
    movetype_sum_stream = new ofstream {m_params.m_output_filebase + ".moves"};
    for (auto& movetype: m_movetypes) {
        movetype->write_log_summary_header(m_logging_stream);
        movetype->write_log_summary(movetype_sum_stream);
        *m_logging_stream << "\n";
    }
    movetype_sum_stream->close();
    delete movetype_sum_stream;
}

void GCMCSimulation::setup_vmd_pipe() {
    string output_filebase {m_params.m_output_filebase + "_vmd"};
    vmd_struct_file = new OrigamiVSFOutputFile {
            output_filebase + ".vsf",
            0,
            m_params.m_max_total_staples,
            m_params.m_max_staple_size,
            m_origami_system};

    vmd_coors_file = new OrigamiVCFOutputFile {
            output_filebase + ".vcf",
            0,
            m_params.m_max_total_staples,
            m_params.m_max_staple_size,
            m_origami_system};

    vmd_states_file = new OrigamiStateOutputFile {
            output_filebase + ".states",
            0,
            m_params.m_max_total_staples,
            m_params.m_max_staple_size,
            m_origami_system};

    vmd_ores_file = new OrigamiOrientationOutputFile {
            output_filebase + ".ores",
            0,
            m_params.m_max_total_staples,
            m_params.m_max_staple_size,
            m_origami_system};

    pipe_to_vmd();
    if (m_params.m_create_vmd_instance) {
        vmd_proc = new bp::child {
                bp::search_path("vmd"),
                "-e",
                m_params.m_vmd_file_dir + "/pipe.tcl",
                "-args",
                m_params.m_vmd_file_dir,
                m_params.m_output_filebase + "_vmd",
                bp::std_out > "/dev/null"};
    }
}

void GCMCSimulation::pipe_to_vmd() {
    vmd_struct_file->open_write_close();
    vmd_coors_file->open_write_close();
    vmd_states_file->open_write_close();
    vmd_ores_file->open_write_close();
}

void GCMCSimulation::close_vmd_pipe() {
    delete vmd_struct_file;
    delete vmd_coors_file;
    delete vmd_states_file;
    delete vmd_ores_file;
}

void GCMCSimulation::close_output_files() {
    for (auto output_file: m_output_files) {
        delete output_file;
    }
    for (auto file: m_config_per_move_files) {
        delete file;
    }
    m_output_files.clear();
}
} // namespace simulation
