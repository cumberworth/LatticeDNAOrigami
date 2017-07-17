// simulation.cpp

#include <random>
#include <iostream>
#include <set>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <sstream>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/process.hpp>

#include "utility.h"
#include "random_gens.h"
#include "movetypes.h"
#include "orientation_movetype.h"
#include "met_movetypes.h"
#include "cb_movetypes.h"
#include "simulation.h"
#include "files.h"

#include "random_gens.h"

namespace Simulation {

    using std::cout;
    using std::min;
    using std::pow;

    namespace mpi = boost::mpi;
    namespace bp = boost::process;

    using namespace Movetypes;
    using namespace OrientationMovetype;
    using namespace MetMovetypes;
    using namespace CBMovetypes;
    using namespace Utility;
    using namespace RandomGen;
    using namespace Files;
    using namespace RandomGen;

    vector<OrigamiOutputFile*> setup_output_files(
            InputParameters& params, string output_filebase,
            OrigamiSystem& origami) {

        // Hack to get a vsf file
        OrigamiVSFOutputFile vsf_file {
                output_filebase + ".vsf", 0,
                params.m_max_total_staples, origami};
        vsf_file.write(0);

        vector<OrigamiOutputFile*> outs {};
        if (params.m_configs_output_freq != 0) {
            OrigamiOutputFile* config_out = new OrigamiTrajOutputFile {
                    output_filebase + ".trj", params.m_configs_output_freq,
                    params.m_max_total_staples, origami};
            outs.push_back(config_out);
        }
        if (params.m_vtf_output_freq != 0) {
            OrigamiOutputFile* config_out = new OrigamiVCFOutputFile {
                    output_filebase + ".vcf", params.m_vtf_output_freq,
                    params.m_max_total_staples, origami};
            outs.push_back(config_out);

            config_out = new OrigamiStateOutputFile {
                    output_filebase + ".states", params.m_vtf_output_freq,
                    params.m_max_total_staples, origami};
            outs.push_back(config_out);

            config_out = new OrigamiOrientationOutputFile {
                    output_filebase + ".ores", params.m_vtf_output_freq,
                    params.m_max_total_staples, origami};
            outs.push_back(config_out);
        }
        if (params.m_counts_output_freq != 0) {
            OrigamiOutputFile* counts_out = new OrigamiCountsOutputFile {
                    output_filebase + ".counts", params.m_counts_output_freq,
                    params.m_max_total_staples, origami};
            outs.push_back(counts_out);
        }
        if (params.m_energies_output_freq != 0) {
            OrigamiOutputFile* energies_out = new OrigamiEnergiesOutputFile {
                    output_filebase + ".ene", params.m_energies_output_freq,
                    params.m_max_total_staples, origami};
            outs.push_back(energies_out);
        }
        if (params.m_order_params_output_freq != 0) {
            OrigamiOutputFile* order_params_out = new OrigamiOrderParamsOutputFile {
                    output_filebase + ".order_params",
                    params.m_order_params_output_freq, params.m_max_total_staples,
                    origami};
            outs.push_back(order_params_out);
        }

        return outs;
    }

    GCMCSimulation::GCMCSimulation(OrigamiSystem& origami_system,
            InputParameters& params) :
            m_origami_system {origami_system},
            m_params {params} {

        m_logging_freq = params.m_logging_freq;
        m_centering_freq = params.m_centering_freq;
        m_constraint_check_freq = params.m_constraint_check_freq;
        m_vmd_pipe_freq = params.m_vmd_pipe_freq;

        if (m_vmd_pipe_freq != 0) {
            setup_vmd_pipe();
        }

        // Constructor movetypes
        construct_movetypes(params);

        // Create cumulative probability array
        double cum_prob {0};
        for (size_t i {0}; i != m_movetypes.size(); i++) {
            cum_prob += params.m_movetype_probs[i];
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
        for (auto movetype_id: params.m_movetypes) {
            shared_ptr<MCMovetype> movetype;
            if (movetype_id == MovetypeID::OrientationRotation) {
                movetype.reset(new OrientationRotationMCMovetype {
                        m_origami_system, m_random_gens, m_ideal_random_walks, params});
            }
            else if (movetype_id == MovetypeID::MetStapleExchange) {
                movetype.reset(new MetStapleExchangeMCMovetype {
                        m_origami_system, m_random_gens, m_ideal_random_walks, params});
            }
            else if (movetype_id == MovetypeID::MetStapleRegrowth) {
                movetype.reset(new MetStapleRegrowthMCMovetype {
                        m_origami_system, m_random_gens, m_ideal_random_walks, params});
            }
            else if (movetype_id == MovetypeID::CBStapleRegrowth) {
                movetype.reset(new CBStapleRegrowthMCMovetype {
                        m_origami_system, m_random_gens, m_ideal_random_walks, params});
            }
            else if (movetype_id == MovetypeID::CTCBScaffoldRegrowth) {
                movetype.reset(new CTCBScaffoldRegrowthMCMovetype {
                        m_origami_system, m_random_gens, m_ideal_random_walks, params});
            }
            else if (movetype_id == MovetypeID::CTCBLinkerRegrowth) {
                movetype.reset(new CTCBLinkerRegrowthMCMovetype {
                        m_origami_system, m_random_gens, m_ideal_random_walks, params});
            }
            m_movetypes.push_back(movetype);
        }
    }

    void GCMCSimulation::simulate(long long int steps, long long int start_step) {

        for (long long int step {start_step + 1}; step != (steps + start_step + 1); step ++) {
            
            // Pick movetype and apply
            shared_ptr<MCMovetype> movetype {select_movetype()};
            bool accepted;
            accepted = movetype->attempt_move();
            if (not accepted) {
                movetype->reset_origami();
            }
            movetype->reset_internal();

            // Center and check constraints
            if (m_centering_freq != 0 and step % m_centering_freq == 0) {
                m_origami_system.centre();
            }
            if (m_constraint_check_freq != 0 and step % m_constraint_check_freq == 0) {
                m_origami_system.check_all_constraints();
            }

            // Write log entry to standard out
            if (m_logging_freq !=0 and step % m_logging_freq == 0) {
                write_log_entry(step, *movetype, accepted);
            }

            // VMD pipe
            if (m_vmd_pipe_freq != 0 and step % m_vmd_pipe_freq == 0) {
                pipe_to_vmd();
            }

            // Update internal simulation variables
            update_internal(step);

            // Write system properties to file
            for (auto output_file: m_output_files) {
                if (output_file->m_write_freq != 0 and step % output_file->m_write_freq == 0) {
                    output_file->write(step);
                }
            }
        }
    }

    shared_ptr<MCMovetype> GCMCSimulation::select_movetype() {
        shared_ptr<MCMovetype> movetype;
        double prob {m_random_gens.uniform_real()};
        for (size_t i {0}; i != m_cumulative_probs.size(); i++) {
            if (prob < m_cumulative_probs[i]) {
                movetype = m_movetypes[i];
                break;
            }
        }
        return movetype;
    }

    void GCMCSimulation::write_log_entry(long long int step, MCMovetype& movetype,
            bool accepted) {

        *m_logging_stream << "Step: " << step << " ";
        *m_logging_stream << "Movetype: " << movetype.m_label() << " ";
        *m_logging_stream << "Staples: " << m_origami_system.num_staples() << " ";
        *m_logging_stream << "Accepted: " << accepted << " ";
        *m_logging_stream << "Temp: " << m_origami_system.m_temp << " ";
        *m_logging_stream << "Energy: " << m_origami_system.energy() << " ";
        *m_logging_stream << "Bias: " << m_origami_system.bias() << " ";
        *m_logging_stream << "\n";
    }

    void GCMCSimulation::setup_vmd_pipe() {
        string output_filebase {m_params.m_output_filebase + "_vmd"};
        vmd_struct_file = new OrigamiVSFOutputFile {
                output_filebase + ".vsf", 0,
                m_params.m_max_total_staples, m_origami_system};

        vmd_coors_file = new OrigamiVCFOutputFile {
                output_filebase + ".vcf", 0,
                m_params.m_max_total_staples, m_origami_system};

        vmd_states_file = new OrigamiStateOutputFile {
                output_filebase + ".states", 0,
                m_params.m_max_total_staples, m_origami_system};

        vmd_ores_file = new OrigamiOrientationOutputFile {
                output_filebase + ".ores", 0,
                m_params.m_max_total_staples, m_origami_system};

        pipe_to_vmd();
        if (m_params.m_create_vmd_instance) {
            vmd_proc = new bp::child {bp::search_path("vmd"),
                "-e", m_params.m_vmd_file_dir + "/pipe.tcl",
                "-args", m_params.m_vmd_file_dir, m_params.m_output_filebase + "_vmd",
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
        m_output_files.clear();
    }
}
