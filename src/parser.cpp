// parser.cpp

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <sstream>
#include <sstream>

#include "parser.h"


namespace Parser {

    namespace po = boost::program_options;
    using std::ifstream;
    using std::cout;

    // Couldn't easily template cause function changes (stoi vs stod)
    vector<int> string_to_int_vector(string string_v) {
        string delim {" "};
        size_t start_pos {0};
        size_t delim_pos {string_v.find(delim)};
        vector<int> v {};
        while (delim_pos != string::npos) {
            v.push_back(stoi(string_v.substr(start_pos, delim_pos)));
            start_pos = delim_pos + 1;
            delim_pos = string_v.find(delim, start_pos);
        }
        v.push_back(stod(string_v.substr(start_pos, string_v.size())));

        return v;
    }

    vector<double> string_to_double_vector(string string_v) {
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

    vector<string> string_to_string_vector(string string_v) {
        std::stringstream sstream {string_v};
        vector<string> v {};
        while (not sstream.eof()) {
            string ele;
            sstream >> ele;
            v.push_back(ele);
        }

        return v;
    }

    template<typename Out>
    void split(const string &s, char delim, Out result) {
        std::stringstream ss;
        ss.str(s);
        string item;
        while (std::getline(ss, item, delim)) {
            *(result++) = item;
        }
    }

    vector<string> split(const string &s, char delim) {
        vector<string> elems;
        Parser::split(s, delim, std::back_inserter(elems));
        return elems;
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

        // Displayed options
        po::options_description displayed_options {"Allowed options"};

        // Command line options description
        po::options_description cl_options {"Command line options"};
        cl_options.add_options()
            ("help,h",
                "Display available options")
            ("parameter_filename,i",
                po::value<string>(),
                "Input file")
        ;
        displayed_options.add(cl_options);

        // Parameter file options description
        po::options_description inp_options {"System input and parameters"};
        inp_options.add_options()
            ("origami_input_filename",
                po::value<string>(&m_origami_input_filename)->default_value(""),
                "Origami input filename")
            ("temp",
                po::value<double>(&m_temp)->default_value(300),
                "System temperature (K)")
            ("staple_M",
                po::value<double>(&m_staple_M)->default_value(1),
                "Staple concentration (mol/L)")
            ("cation_M",
                po::value<double>(&m_cation_M)->default_value(1),
                "Cation concentration (mol/L)")
            ("temp_for_staple_u",
                po::value<double>(&m_temp_for_staple_u)->default_value(300),
                "Temperature to calculate chemical potential with")
            ("staple_u_mult",
                po::value<double>(&m_staple_u_mult)->default_value(1),
                "Multiplier for staple u")
            ("lattice_site_volume",
                po::value<double>(&m_lattice_site_volume)->default_value(1),
                "Volume per lattice site (L)")
            ("cyclic",
                po::value<bool>(&m_cyclic)->default_value(false),
                "Cyclic scaffold")
            ("no_misbinding",
                po::value<bool>(&m_no_misbinding)->default_value(false),
                "Turn off misbinding")
            ("energy_filebase",
                po::value<string>(&m_energy_filebase)->default_value(""),
                "Filebase for read/write of energies")
            ("restart_traj_file",
                po::value<string>(&m_restart_traj_file)->default_value(""),
                "Trajectory file to restart from")
            ("restart_traj_filebase",
                po::value<string>(&m_restart_traj_filebase)->default_value(""),
                "Trajectory restart filebase")
            ("restart_traj_postfix",
                po::value<string>(&m_restart_traj_postfix)->default_value(".trj"),
                "Trajectory restart postfix")
            ("restart_traj_files",
                po::value<string>(),
                "Trajectory restart files for each replicate")
            ("restart_step",
                po::value<int>(&m_restart_step)->default_value(0),
                "Step to restart from")
            ("restart_steps",
                po::value<string>(),
                "Restart step for each replicate")
            ("vmd_file_dir",
                po::value<string>(&m_vmd_file_dir)->default_value(""),
                "Directory containing VMD scripts for viewing simulations")
            ;
        displayed_options.add(inp_options);

        po::options_description op_options {"Order parameter options"};
        op_options.add_options()
            ("distance_pairs",
                po::value<string>(),
                "Scaffold domains to restrain")
            ("distance_sum",
                po::value<bool>(&m_distance_sum)->default_value(false),
                "Sum of distance pairs")
            ;
        displayed_options.add(op_options);

        po::options_description bias_options {"Bias function options"};
        bias_options.add_options()
            ("distance_bias",
                po::value<bool>(&m_distance_bias)->default_value(false),
                "Include domain pair distance bias")
            ("min_dist",
                po::value<int>(&m_min_dist)->default_value(0),
                "Distance at which bias switched off")
            ("max_dist",
                po::value<int>(&m_max_dist)->default_value(0),
                "Distance at which bias switched on")
            ("max_bias",
                po::value<double>(&m_max_bias)->default_value(0),
                "Value of bias at dist >= max_dist")
            ("bias_mult",
                po::value<double>(&m_bias_mult)->default_value(1),
                "Total bias multiplier")
            ("grid_bias",
                po::value<bool>(&m_grid_bias)->default_value(false),
                "Include a grid bias")
            ("square_well_bias",
                po::value<bool>(&m_square_well_bias)->default_value(false),
                "Include a square well bias on the number of staples")
            ("min_well_param",
                po::value<int>(&m_min_well_param)->default_value(0),
                "Min number of staples for well bias")
            ("max_well_param",
                po::value<int>(&m_max_well_param)->default_value(0),
                "Max number of staples for well bias ")
            ("well_bias",
                po::value<double>(&m_well_bias)->default_value(0),
                "Well bias")
            ("outside_bias",
                po::value<double>(&m_outside_bias)->default_value(0),
                "Bias outside the well")
        ;
        displayed_options.add(bias_options);

        po::options_description sim_options {"General simulation options"};
        sim_options.add_options()
            ("simulation_type",
                po::value<string>(&m_simulation_type)->default_value("constant_temp"),
                "constant_temp, annealing, or parallel_tempering")
            ("steps",
                po::value<long long int>(&m_steps)->default_value(0),
                "Number of MC steps")
            ("logging_freq",
                po::value<int>(&m_logging_freq)->default_value(0),
                "Logging frequency")
            ("centering_freq",
                po::value<int>(&m_centering_freq)->default_value(0),
                "Centering frequency")
            ("centering_domain",
                po::value<int>(&m_centering_domain)->default_value(0),
                "Domain to center on")
            ("constraint_check_freq",
                po::value<int>(&m_constraint_check_freq)->default_value(0),
                "Constraint check frequency")
            ("orientation_rotation",
                po::value<string>(),
                "Orientational rotation movetype probability")
            ("met_staple_exchange",
                po::value<string>(),
                "Met staple exchange movetype probability")
            ("met_staple_regrowth",
                po::value<string>(),
                "Met staple regrowth movetype probability")
            ("cb_staple_exchange",
                po::value<string>(),
                "CB staple exchange movetype probability")
            ("cb_staple_regrowth",
                po::value<string>(),
                "CB staple regrowth movetype probability")
            ("ctcb_scaffold_regrowth",
                po::value<string>(),
                "CTCB scaffold regrowth movetype probability")
            ("ctcb_linker_regrowth",
                po::value<string>(),
                "CTCB linker regrowth movetype probability")
            ("max_displacement",
                po::value<int>(&m_max_displacement)->default_value(0),
                "Max displacement of linked region in CTCB linker regrowth")
            ("max_turns",
                po::value<int>(&m_max_turns)->default_value(0),
                "Max number of quarter turns of linked region in CTCB linker regrowth")
            ("num_walks_filename",
                po::value<string>(&m_num_walks_filename)->default_value(""),
                "Precalculated number of ideal random walks archive")
            ("exchange_mult",
                po::value<int>(&m_exchange_mult)->default_value(1),
                "Exchange acceptance probability multiplier")
            ("max_total_staples",
                po::value<int>(&m_max_total_staples)->default_value(999),
                "Max number of total staples")
            ("max_type_staples",
                po::value<int>(&m_max_type_staples)->default_value(999),
                "Max number of staples of a given type")
        ;
        displayed_options.add(sim_options);

        po::options_description enum_options {"Enumeration options"};
        enum_options.add_options()
            ("enumerate_staples_only",
                po::value<bool>(&m_enumerate_staples_only)->default_value(false),
                "Enumerate staples only")
        ;
        displayed_options.add(enum_options);

        po::options_description annealing_options {"Annealing simulation options"};
        annealing_options.add_options()
            ("max_temp",
                po::value<double>(&m_max_temp)->default_value(400),
                "Maximum temperature for annealing")
            ("min_temp",
                po::value<double>(&m_min_temp)->default_value(300),
                "Minimum temperature for annealing")
            ("temp_interval",
                po::value<double>(&m_temp_interval)->default_value(1),
                "Temperature interval for annealing")
            ("steps_per_temp",
                po::value<int>(&m_steps_per_temp)->default_value(0),
                "Steps per temperature in annealing")
        ;
        displayed_options.add(annealing_options);

        po::options_description pt_options {"Parallel tempering simulation options"};
        pt_options.add_options()
            ("temps",
                po::value<string>(),
                "Temperature list")
            ("num_reps",
                po::value<int>(&m_num_reps)->default_value(1),
                "Number of replicas")
            ("exchange_interval",
                po::value<int>(&m_exchange_interval)->default_value(0),
                "Steps between exchange attempts")
            ("constant_staple_M",
                po::value<bool>(&m_constant_staple_M)->default_value(true),
                "Hold staple concentration constant")
            ("chem_pot_mults",
                po::value<string>(),
                "Factor to multiply base chem pot for each rep")
            ("bias_mults",
                po::value<string>(),
                "Multiplier for system bias")
            ("restart_swap_file",
                po::value<string>(&m_restart_swap_file)->default_value(""),
                "Swap file to restart from")
        ;
        displayed_options.add(pt_options);

        po::options_description us_options {"Umbrella sampling simulation options"};
        us_options.add_options()
            ("order_params",
                po::value<string>(),
                "List of order parameters to apply to")
            ("max_num_iters",
                po::value<int>(&m_max_num_iters)->default_value(0),
                "Number of iterations")
            ("max_D_bias",
                po::value<double>(&m_max_D_bias)->default_value(0),
                "Max change in bias per iteration")
            ("equil_steps",
                po::value<long long int>(&m_equil_steps)->default_value(0),
                "Number of equilibration steps")
            ("iter_steps",
                po::value<long long int>(&m_iter_steps)->default_value(0),
                "Number of steps per iteration")
            ("prod_steps",
                po::value<long long int>(&m_prod_steps)->default_value(0),
                "Number of production steps")
            ("max_rel_P_diff",
                po::value<double>(&m_max_rel_P_diff)->default_value(0.1),
                "Maximum allowed change in P for convergence")
            ("biases_file",
                po::value<string>(&m_biases_file)->default_value(""),
                "Initial guesses at grid biases")
            ("biases_filebase",
                po::value<string>(&m_biases_filebase)->default_value(""),
                "Filebase for grid bias files")
            ("multi_window",
                po::value<bool>(&m_multi_window)->default_value(false),
                "Use multiple windows")
            ("windows_file",
                po::value<string>(&m_windows_file)->default_value(""),
                "File containing windows as min/max pairs of tuples")
        ;
        displayed_options.add(us_options);

        po::options_description out_options {"Output options"};
        out_options.add_options()
            ("output_filebase",
                po::value<string>(&m_output_filebase)->default_value(""),
                "Base name for output files")
            ("configs_output_freq",
                po::value<int>(&m_configs_output_freq)->default_value(0),
                "Configuration output write frequency")
            ("vtf_output_freq",
                po::value<int>(&m_vtf_output_freq)->default_value(0),
                "Configuration output write frequency")
            ("vcf_per_domain",
                po::value<bool>(&m_vcf_per_domain)->default_value(false),
                "Write a VCF entry for every domain grown")
            ("counts_output_freq",
                po::value<int>(&m_counts_output_freq)->default_value(0),
                "Counts output write frequency")
            ("energies_output_freq",
                po::value<int>(&m_energies_output_freq)->default_value(0),
                "Energies output write frequency")
            ("order_params_output_freq",
                po::value<int>(&m_order_params_output_freq)->default_value(0),
                "Order parameters write frequency")
            ("vmd_pipe_freq",
                po::value<int>(&m_vmd_pipe_freq),
                "Realtime VMD visualization updating frequency")
            ("create_vmd_instance",
                po::value<bool>(&m_create_vmd_instance)->default_value(false),
                "Create VMD instance")
        ;
        displayed_options.add(out_options);

        // Command line variable map
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, displayed_options), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << "\n";
            cout << displayed_options;
            cout << "\n";
            exit(1);
        }
        if (not vm.count("parameter_filename")) {
            cout << "Input parameter file must be provided.\n";
            exit(1);
        }

        // Parameter file variable map
        ifstream parameter_file {vm["parameter_filename"].as<string>()};
        po::store(po::parse_config_file(parameter_file, displayed_options), vm);
        po::notify(vm);

        set_default_sim_options();

        // Process custom vector and custom types
        if (vm.count("restart_traj_files")) {
            string file_s = vm["restart_traj_files"].as<string>();
            m_restart_traj_files = string_to_string_vector(file_s);
        }
        if (vm.count("distance_pairs")) {
            string distance_pairs_s= vm["distance_pairs"].as<string>();
            m_distance_pairs = string_to_int_vector(distance_pairs_s);
        }
        if (vm.count("distance_bias")) {
            m_distance_bias = vm["distance_bias"].as<bool>();
            m_biases_present = true;
        }
        if (vm.count("max_bias")) {
            m_max_bias = vm["max_bias"].as<double>();
            m_biases_present = true;
        }
        if (vm.count("grid_bias")) {
            m_grid_bias = vm["grid_bias"].as<bool>();
            m_biases_present = true;
        }
        if (vm.count("square_well_bias")) {
            m_square_well_bias = vm["square_well_bias"].as<bool>();
            m_biases_present = true;
        }
        if (vm.count("orientation_rotation")) {
            string unparsed_fraction {vm["orientation_rotation"].as<string>()};
            Fraction prob {unparsed_fraction};
            m_movetype_probs.push_back(prob.to_double());
            m_movetypes.push_back(MovetypeID::OrientationRotation);
        }
        if (vm.count("met_staple_exchange")) {
            string unparsed_fraction {vm["met_staple_exchange"].as<string>()};
            Fraction prob {unparsed_fraction};
            m_movetype_probs.push_back(prob.to_double());
            m_movetypes.push_back(MovetypeID::MetStapleExchange);
        }
        if (vm.count("met_staple_regrowth")) {
            string unparsed_fraction {vm["met_staple_regrowth"].as<string>()};
            Fraction prob {unparsed_fraction};
            m_movetype_probs.push_back(prob.to_double());
            m_movetypes.push_back(MovetypeID::MetStapleRegrowth);
        }
        if (vm.count("cb_staple_regrowth")) {
            string unparsed_fraction {vm["cb_staple_regrowth"].as<string>()};
            Fraction prob {unparsed_fraction};
            m_movetype_probs.push_back(prob.to_double());
            m_movetypes.push_back(MovetypeID::CBStapleRegrowth);
        }
        if (vm.count("ctcb_scaffold_regrowth")) {
            string unparsed_fraction {vm["ctcb_scaffold_regrowth"].as<string>()};
            Fraction prob {unparsed_fraction};
            m_movetype_probs.push_back(prob.to_double());
            m_movetypes.push_back(MovetypeID::CTCBScaffoldRegrowth);
        }
        if (vm.count("ctcb_linker_regrowth")) {
            string unparsed_fraction {vm["ctcb_linker_regrowth"].as<string>()};
            Fraction prob {unparsed_fraction};
            m_movetype_probs.push_back(prob.to_double());
            m_movetypes.push_back(MovetypeID::CTCBLinkerRegrowth);
        }
        if (vm.count("temps")) {
            string temps_s = vm["temps"].as<string>();
            m_temps = string_to_double_vector(temps_s);
        }
        if (vm.count("chem_pot_mults")) {
            string chem_pot_mults_s = vm["chem_pot_mults"].as<string>();
            m_chem_pot_mults = string_to_double_vector(chem_pot_mults_s);
        }
        if (vm.count("bias_mults")) {
            string bias_mults_s {vm["bias_mults"].as<string>()};
            m_bias_mults = string_to_double_vector(bias_mults_s);
        }
        if (vm.count("order_params")) {
            string order_params_s {vm["order_params"].as<string>()};
            m_order_params = split(order_params_s, ' ');
        }
    }

    void InputParameters::set_default_sim_options() {
        if (m_simulation_type == "umbrella_sampling" or
                m_simulation_type == "mw_umbrella_sampling") {
            m_biases_present = true;
            m_grid_bias = true;
        }
    }
}
