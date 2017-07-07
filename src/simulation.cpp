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
        // Setup trajectory and staple/domain count file

        vector<OrigamiOutputFile*> outs {};
        if (params.m_configs_output_freq != 0) {
            OrigamiOutputFile* config_out = new OrigamiTrajOutputFile {
                    output_filebase + ".trj", params.m_configs_output_freq,
                    origami};
            outs.push_back(config_out);
        }
        if (params.m_counts_output_freq != 0) {
            OrigamiOutputFile* counts_out = new OrigamiCountsOutputFile {
                    output_filebase + ".counts", params.m_counts_output_freq,
                    origami};
            outs.push_back(counts_out);
        }
        if (params.m_energies_output_freq != 0) {
            OrigamiOutputFile* energies_out = new OrigamiEnergiesOutputFile {
                    output_filebase + ".ene", params.m_energies_output_freq,
                    origami};
            outs.push_back(energies_out);
        }
        if (params.m_order_params_output_freq != 0) {
            OrigamiOutputFile* order_params_out = new OrigamiOrderParamsOutputFile {
                    output_filebase + ".order_params",
                    params.m_order_params_output_freq, origami};
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
                m_origami_system.check_all_constraints();
            }

            // Write log entry to standard out
            if (m_logging_freq !=0 and step % m_logging_freq == 0) {
                write_log_entry(step, *movetype, accepted);
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

    void GCMCSimulation::close_output_files() {
        for (auto output_file: m_output_files) {
            delete output_file;
        }
        m_output_files.clear();
    }
}
