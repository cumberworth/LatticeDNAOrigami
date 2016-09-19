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
        ("origami_input_filename", po::value<string>(), "Origami input filename")
        ("configs_output_filename", po::value<string>(), "Configuration output filename")
        ("configs_output_freq", po::value<int>(), "Configuration output write frequency")
        ("temp", po::value<double>(), "System temperature (K)")
        ("staple_M", po::value<double>(), "Staple concentration (mol/L)")
        ("cation_M", po::value<double>(), "Cation concentration (mol/L)")
        ("lattice_site_volume", po::value<double>(), "Volume per lattice site (L)")
        ("cyclic", po::value<bool>(), "Cyclic scaffold")
        ("steps", po::value<int>(), "Number of MC steps")
        ("logging_freq", po::value<int>(), "Logging frequency")
        ("centering_freq", po::value<int>(), "Centering frequency")
        ("orientation_rotation", po::value<double>(), "Orientational rotation movetype probability")
        ("cb_staple_exchange", po::value<double>(), "CB staple exchange movetype probability")
        ("cb_staple_regrowth", po::value<double>(), "CB staple regrowth movetype probability")
        ("ctcb_scaffold_regrowth", po::value<double>(), "CTCB scaffold regrowth movetype probability")
        ;

    po::variables_map vm;
    po::store(po::parse_config_file(parameter_file, desc), vm);
    po::notify(vm);

    // Extract variables from variables map
    if (vm.count("origami_input_filename")) {
        m_origami_input_filename = vm["origami_input_filename"].as<string>();
    }
    else {
        cout << "Origami input filename needed.";
        exit(1);
    }

    if (vm.count("configs_output_freq")) {
        m_configs_output_freq = vm["configs_output_freq"].as<int>();
    }

    if (m_configs_output_freq != 0) {
        if (vm.count("configs_output_filename"))
            m_configs_output_filename = vm["configs_output_filename"].as<string>();
        else {
            cout << "Configuration output filename needed.";
            exit(1);
        }
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
        m_movetype_probs.push_back(vm["orientation_rotation"].as<double>());
        m_movetype_constructors.push_back(movetype[1]);
    }

    if (vm.count("cb_staple_exchange")) {
        m_movetype_probs.push_back(vm["cb_staple_exchange"].as<double>());
        m_movetype_constructors.push_back(movetype[3]);
    }

    if (vm.count("cb_staple_regrowth")) {
        m_movetype_probs.push_back(vm["cb_staple_regrowth"].as<double>());
        m_movetype_constructors.push_back(movetype[4]);
    }

    if (vm.count("ctcb_scaffold_regrowth")) {
        m_movetype_probs.push_back(vm["ctcb_scaffold_regrowth"].as<double>());
        m_movetype_constructors.push_back(movetype[5]);
    }
}
