// parser.cpp

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

#include "parser.h"
#include "movetypes.h"

namespace po = boost::program_options;

using std::ifstream;
using std::cout;

using namespace Parser;
using namespace Movetypes;

vector<double> Parser::string_to_double_vector(string string_v) {
    string delim {" "};
    size_t start_pos {0};
    size_t delim_pos {string_v.find(delim)};
    vector<double> v {};
    while (delim_pos != string::npos) {
        v.push_back(stod(string_v.substr(start_pos, delim_pos)));
        start_pos = delim_pos + 1;
        delim_pos = string_v.find(delim, start_pos);
    }
    v.push_back(stod(string_v.substr(start_pos, string_v.size())));

    return v;
}

Fraction::Fraction(string unparsed_fraction) {
    string delimiter {"/"};
    auto delim_pos {unparsed_fraction.find(delimiter)};

    // Assume is real number
    if (delim_pos == string::npos) {
        m_numerator = stod(unparsed_fraction);
        m_denominator = 1;
    }
    else {
        auto end_pos {unparsed_fraction.size()};
        m_numerator = stod(unparsed_fraction.substr(0, delim_pos));
        m_denominator = stod(unparsed_fraction.substr(delim_pos + 1, end_pos));
    }
    m_double_fraction = m_numerator / m_denominator;
}

InputParameters::InputParameters(int argc, char* argv[]) {

    // Command line parser for getting config file name
    po::options_description cl_desc {"Allowed options"};
    cl_desc.add_options()
        ("parameter_filename,i", po::value<string>())
        ;

    po::variables_map cl_vm;
    po::store(po::parse_command_line(argc, argv, cl_desc), cl_vm);
    po::notify(cl_vm);

    if (not cl_vm.count("parameter_filename")) {
        cout << "Input parameter file must be provided.\n";
        exit(1);
    }
    ifstream parameter_file {cl_vm["parameter_filename"].as<string>()};

    // Setup parser
    po::options_description desc {};
    desc.add_options()

        // System input and parameters
        ("origami_input_filename", po::value<string>(), "Origami input filename")
        ("temp", po::value<double>(), "System temperature (K)")
        ("staple_M", po::value<double>(), "Staple concentration (mol/L)")
        ("cation_M", po::value<double>(), "Cation concentration (mol/L)")
        ("lattice_site_volume", po::value<double>(), "Volume per lattice site (L)")
        ("cyclic", po::value<bool>(), "Cyclic scaffold")

        // General simulation parameters
        ("simulation_type", po::value<string>(), "constant_temp, annealing, or parallel_tempering")
        ("steps", po::value<int>(), "Number of MC steps")
        ("logging_freq", po::value<int>(), "Logging frequency")
        ("centering_freq", po::value<int>(), "Centering frequency")
        ("orientation_rotation", po::value<string>(), "Orientational rotation movetype probability")
        ("met_staple_exchange", po::value<string>(), "Met staple exchange movetype probability")
        ("met_staple_regrowth", po::value<string>(), "Met staple regrowth movetype probability")
        ("cb_staple_exchange", po::value<string>(), "CB staple exchange movetype probability")
        ("cb_staple_regrowth", po::value<string>(), "CB staple regrowth movetype probability")
        ("ctcb_scaffold_regrowth", po::value<string>(), "CTCB scaffold regrowth movetype probability")

        // Annealing simulation parameters
        ("max_temp", po::value<double>(), "Maximum temperature for annealing")
        ("min_temp", po::value<double>(), "Minimum temperature for annealing")
        ("temp_interval", po::value<double>(), "Temperature interval for annealing")
        ("steps_per_temp", po::value<int>(), "Steps per temperature in annealing")

        // Parallel tempering options 
        ("temps", po::value<string>(), "Temperature list")
        ("num_reps", po::value<int>(), "Number of replicas")
        ("exchange_interval", po::value<int>(), "Steps between exchange attempts")


        // Output options
        ("output_filebase", po::value<string>(), "Base name for output files")
        ("configs_output_freq", po::value<int>(), "Configuration output write frequency")
        ("counts_output_freq", po::value<int>(), "Coounts output write frequency")
        ;

    po::variables_map vm;
    po::store(po::parse_config_file(parameter_file, desc), vm);
    po::notify(vm);

    // Extract variables from variables map
    if (vm.count("origami_input_filename")) {
        m_origami_input_filename = vm["origami_input_filename"].as<string>();
    }

    if (vm.count("configs_output_freq")) {
        m_configs_output_freq = vm["configs_output_freq"].as<int>();
    }

    if (vm.count("counts_output_freq")) {
        m_counts_output_freq = vm["counts_output_freq"].as<int>();
    }

    if (vm.count("temp")) {
        m_temp = vm["temp"].as<double>();
    }

    if (vm.count("staple_M")) {
        m_staple_M = vm["staple_M"].as<double>();
    }

    if (vm.count("cation_M")) {
        m_cation_M = vm["cation_M"].as<double>();
    }

    if (vm.count("lattice_site_volume")) {
        m_lattice_site_volume = vm["lattice_site_volume"].as<double>();
    }

    if (vm.count("cyclic")) {
        m_cyclic = vm["cyclic"].as<bool>();
    }

    if (vm.count("steps")) {
        m_steps = vm["steps"].as<int>();
    }

    if (vm.count("logging_freq")) {
        m_logging_freq = vm["logging_freq"].as<int>();
    }

    if (vm.count("centering_freq")) {
        m_centering_freq = vm["centering_freq"].as<int>();
    }

    if (vm.count("orientation_rotation")) {
        string unparsed_fraction {vm["orientation_rotation"].as<string>()};
        Fraction prob {unparsed_fraction};
        m_movetype_probs.push_back(prob.to_double());
        m_movetype_constructors.push_back(movetype[1]);
    }

    if (vm.count("met_staple_exchange")) {
        string unparsed_fraction {vm["met_staple_exchange"].as<string>()};
        Fraction prob {unparsed_fraction};
        m_movetype_probs.push_back(prob.to_double());
        m_movetype_constructors.push_back(movetype[2]);
    }

    if (vm.count("met_staple_regrowth")) {
        string unparsed_fraction {vm["met_staple_regrowth"].as<string>()};
        Fraction prob {unparsed_fraction};
        m_movetype_probs.push_back(prob.to_double());
        m_movetype_constructors.push_back(movetype[3]);
    }

    if (vm.count("cb_staple_exchange")) {
        string unparsed_fraction {vm["cb_staple_exchange"].as<string>()};
        Fraction prob {unparsed_fraction};
        m_movetype_probs.push_back(prob.to_double());
        m_movetype_constructors.push_back(movetype[4]);
    }

    if (vm.count("cb_staple_regrowth")) {
        string unparsed_fraction {vm["cb_staple_regrowth"].as<string>()};
        Fraction prob {unparsed_fraction};
        m_movetype_probs.push_back(prob.to_double());
        m_movetype_constructors.push_back(movetype[5]);
    }

    if (vm.count("ctcb_scaffold_regrowth")) {
        string unparsed_fraction {vm["ctcb_scaffold_regrowth"].as<string>()};
        Fraction prob {unparsed_fraction};
        m_movetype_probs.push_back(prob.to_double());
        m_movetype_constructors.push_back(movetype[6]);
    }

    if (vm.count("simulation_type")) {
        m_simulation_type = vm["simulation_type"].as<string>();
    }

    if (vm.count("max_temp")) {
        m_max_temp = vm["max_temp"].as<double>();
    }

    if (vm.count("min_temp")) {
        m_min_temp = vm["min_temp"].as<double>();
    }

    if (vm.count("temp_interval")) {
        m_temp_interval = vm["temp_interval"].as<double>();
    }

    if (vm.count("steps_per_temp")) {
        m_steps_per_temp = vm["steps_per_temp"].as<int>();
    }

    if (vm.count("temps")) {
        string temps_s = vm["temps"].as<string>();
        m_temps = string_to_double_vector(temps_s);
    }

    if (vm.count("num_reps")) {
        m_num_reps = vm["num_reps"].as<int>();
    }

    if (vm.count("exchange_interval")) {
        m_exchange_interval = vm["exchange_interval"].as<int>();
    }

    if (vm.count("output_filebase")) {
        m_output_filebase = vm["output_filebase"].as<string>();
    }
}
