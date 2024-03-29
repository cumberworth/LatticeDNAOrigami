// parser.cpp

#include <fstream>
#include <iostream>

#include "LatticeDNAOrigami/parser.hpp"
#include "LatticeDNAOrigami/utility.hpp"
#include "LatticeDNAOrigami/version.hpp"

namespace parser {

using std::cout;
using std::ifstream;

using utility::string_to_double_vector;
using utility::string_to_int_vector;
using utility::string_to_string_vector;

InputParameters::InputParameters(int argc, char* argv[]) {

    // Displayed options
    po::options_description displayed_options {"Allowed options"};

    // Command line options description
    po::options_description cl_options {"Command line options"};
    cl_options.add_options()("help,h", "Display available options")(
            "parameter_filename,i", po::value<string>(), "Input file")(
            "version,v", "Print program version");
    displayed_options.add(cl_options);

    // Parameter file options description
    po::options_description inp_options {"System input and parameters"};
    inp_options.add_options()(

            "origami_input_filename",
            po::value<string>(&m_origami_input_filename),
            "Origami input filename, specifying scaffold and staples")(

            "domain_type",
            po::value<string>(&m_domain_type)->default_value("Halfturn"),
            "Domain type [HalfTurn, ThreeQuarterTurn]")(

            "binding_pot",
            po::value<string>(&m_binding_pot)->default_value("FourBody"),
            "Binding potential to use; only 'FourBody' is defined")(

            "misbinding_pot",
            po::value<string>(&m_misbinding_pot)->default_value("Opposing"),
            "Misbinding potential to use [Opposing, Disallowed]")(

            "stacking_pot",
            po::value<string>(&m_stacking_pot)->default_value("Constant"),
            "Stacking potential to use [Constant, SequenceSpecific]")(

            "hybridization_pot",
            po::value<string>(&m_hybridization_pot)
                    ->default_value("NearestNeighbour"),
            "Hybridization potential to use [NearestNeighbour, Uniform, "
            "Specified]")(

            "apply_mean_field_cor",
            po::value<bool>(&m_apply_mean_field_cor)->default_value(false),
            "Apply mean field correction")(

            "temp",
            po::value<double>(&m_temp)->default_value(300),
            "System temperature (K)")(

            "staple_M",
            po::value<double>(&m_staple_M)->default_value(1),
            "Staple concentration (mol/L)")(

            "cation_M",
            po::value<double>(&m_cation_M)->default_value(1),
            "Cation concentration (mol/L)")(

            "staple_u_mult",
            po::value<double>(&m_staple_u_mult)->default_value(1),
            "Multiplier for staple u")(

            "stacking_ene",
            po::value<double>(&m_stacking_ene)->default_value(1),
            "Stacking energy for constant potential (kb K)")(

            "binding_h",
            po::value<double>(&m_binding_h)->default_value(0),
            "Hybridization enthalpy for bound domains (kb K); used with "
            "'uniform' "
            "hybridzation potential")(

            "binding_s",
            po::value<double>(&m_binding_s)->default_value(0),
            "Hybridization entropy for bound domains (kb); used with 'uniform' "
            "hybridzation potential")(

            "misbinding_h",
            po::value<double>(&m_misbinding_h)->default_value(0),
            "Hybridization enthalpy for misbound domains (kb K); used with "
            "'uniform' "
            "hybridzation potential")(

            "misbinding_s",
            po::value<double>(&m_misbinding_s)->default_value(0),
            "Hybridization entropy for misbound domains (kb); used with "
            "'uniform' "
            "hybridzation potential")(

            "min_total_staples",
            po::value<int>(&m_min_total_staples)->default_value(0),
            "Min number of total staples")(

            "max_total_staples",
            po::value<int>(&m_max_total_staples)->default_value(999),
            "Max number of total staples")(

            "max_type_staples",
            po::value<int>(&m_max_type_staples)->default_value(999),
            "Max number of staples of a given type")(

            "max_staple_size",
            po::value<int>(&m_max_staple_size)->default_value(2),
            "Max number of domains per staple")(

            "excluded_staples",
            po::value<string>(),
            "Staple types to exclude (this is in developement)")(

            "domain_update_biases_present",
            po::value<bool>(&m_domain_update_biases_present)
                    ->default_value(false),
            "Biases present that require an update after each domain has been "
            "set")(

            "order_parameter_file",
            po::value<string>(&m_ops_filename)->default_value(""),
            "Order parameter specification file")(

            "bias_functions_file",
            po::value<string>(&m_bias_funcs_filename)->default_value(""),
            "Bias function specification file")(

            "bias_functions_mult",
            po::value<double>(&m_bias_funcs_mult)->default_value(1),
            "System bias function multiplier")(

            "energy_filebase",
            po::value<string>(&m_energy_filebase)->default_value(""),
            "Filebase for read/write of energies (currently disabled)")(

            "simulation_type",
            po::value<string>(&m_simulation_type)
                    ->default_value("constant_temp"),
            "Simulation type [enumerate, constant_temp, annealing, "
            "t_parallel_tempering, ut_parallel_tempering, "
            "hut_parallel_tempering, st_parallel_tempering, "
            "2d_parallel_tempering, umbrella_sampling, mw_umbrella_sampling, "
            "ptmw_umbrella_sampling]");

    displayed_options.add(inp_options);

    po::options_description enum_options {"Enumeration options"};
    enum_options.add_options()(

            "enumerate_staples_only",
            po::value<bool>(&m_enumerate_staples_only)->default_value(false),
            "Enumerate staples only");

    displayed_options.add(enum_options);

    po::options_description sim_options {"General simulation options"};
    sim_options.add_options()(

            "random_seed",
            po::value<int>(&m_random_seed)->default_value(-1),
            "Seed for random number generator")(

            "movetype_file",
            po::value<string>(&m_movetype_filename),
            "Movetype specificiation file")(

            "read_num_walks",
            po::value<bool>(&m_read_num_walks)->default_value(false),
            "Read precalculated number of ideal random walks archive")(

            "num_walks_filename",
            po::value<string>(&m_num_walks_filename)->default_value(""),
            "Precalculated number of ideal random walks archive")(

            "restart_from_config",
            po::value<bool>(&m_restart_from_config)->default_value(false),
            "Restart simulation from provided configuration(s)")(

            "restart_traj_file",
            po::value<string>(&m_restart_traj_file)->default_value(""),
            "Trajectory file to restart from")(

            "restart_traj_filebase",
            po::value<string>(&m_restart_traj_filebase)->default_value(""),
            "Trajectory restart filebase")(

            "restart_traj_postfix",
            po::value<string>(&m_restart_traj_postfix)->default_value(".trj"),
            "Trajectory restart postfix")(

            "restart_traj_files",
            po::value<string>(),
            "Trajectory restart files for each replicate")(

            "restart_from_swap",
            po::value<bool>(&m_restart_from_swap)->default_value(false),
            "Restart simulation from swap file")(

            "restart_us_iter",
            po::value<bool>(&m_restart_us_iter)->default_value(false),
            "Restart US simulation iteration")(

            "restart_us_filebase",
            po::value<string>(&m_restart_us_filebase)->default_value(""),
            "Trajectory restart filebase")(

            "restart_step",
            po::value<int>(&m_restart_step)->default_value(0),
            "Step to restart from")(

            "restart_steps",
            po::value<string>(),
            "Restart step for each replicate")(

            "read_rand_engine_state",
            po::value<bool>(&m_read_rand_engine_state)->default_value(false),
            "Read random engine state file")(

            "rand_engine_state_file",
            po::value<string>(&m_rand_engine_state_file)->default_value(""),
            "Random engine state file")(

            "vmd_file_dir",
            po::value<string>(&m_vmd_file_dir)->default_value(""),
            "Directory containing VMD scripts for viewing simulations")(

            "centering_freq",
            po::value<int>(&m_centering_freq)->default_value(0),
            "Centering frequency")(

            "centering_domain",
            po::value<int>(&m_centering_domain)->default_value(0),
            "Domain to center on")(

            "constraint_check_freq",
            po::value<int>(&m_constraint_check_freq)->default_value(0),
            "Constraint check frequency")(

            "allow_nonsensical_ps",
            po::value<bool>(&m_allow_nonsensical_ps)->default_value(false),
            "Allow nonsensical exchange probabilities")(

            "max_duration",
            po::value<double>(&m_max_duration)->default_value(10e9),
            "Maximum duration of simulation (s)");

    displayed_options.add(sim_options);

    po::options_description cons_t_options {"Constant temperature options"};
    cons_t_options.add_options()(

            "ct_steps",
            po::value<long long int>(&m_ct_steps)->default_value(0),
            "Number of MC steps");

    displayed_options.add(cons_t_options);

    po::options_description annealing_options {"Annealing simulation options"};
    annealing_options.add_options()(

            "constant_staple_M",
            po::value<bool>(&m_constant_staple_M)->default_value(true),
            "Hold staple concentration constant")(

            "max_temp",
            po::value<double>(&m_max_temp)->default_value(400),
            "Maximum temperature for annealing (K)")(

            "min_temp",
            po::value<double>(&m_min_temp)->default_value(300),
            "Minimum temperature for annealing (K)")(

            "temp_interval",
            po::value<double>(&m_temp_interval)->default_value(1),
            "Temperature interval for annealing (K)")(

            "steps_per_temp",
            po::value<long long int>(&m_steps_per_temp)->default_value(0),
            "Steps per temperature in annealing");

    displayed_options.add(annealing_options);

    po::options_description pt_options {
            "Parallel tempering simulation options"};
    pt_options.add_options()("temps", po::value<string>(), "Temperature list")(

            "num_reps",
            po::value<int>(&m_num_reps)->default_value(1),
            "Number of replicas")(

            "swaps",
            po::value<long long int>(&m_swaps)->default_value(0),
            "Number of swaps")(

            "max_pt_dur",
            po::value<double>(&m_max_pt_dur)->default_value(10e9),
            "Maximum duration (s)")(

            "exchange_interval",
            po::value<int>(&m_exchange_interval)->default_value(0),
            "Steps between exchange attempts")(

            "chem_pot_mults",
            po::value<string>(),
            "Factors to multiply base chem pot for each rep (one for each "
            "rep)")(

            "bias_mults",
            po::value<string>(),
            "Multipliers for system bias (one for each rep)")(

            "stacking_mults",
            po::value<string>(),
            "Multipliers for stacking energy (one for each rep)")(

            "restart_swap_file",
            po::value<string>(&m_restart_swap_file)->default_value(""),
            "Swap file to restart from");

    displayed_options.add(pt_options);

    po::options_description us_options {"Umbrella sampling simulation options"};
    us_options.add_options()(

            "us_grid_bias_tag",
            po::value<string>(&m_us_grid_bias_tag)->default_value(""),
            "Tag of grid bias function to use for US")(

            "max_num_iters",
            po::value<int>(&m_max_num_iters)->default_value(0),
            "Number of iterations")(

            "max_D_bias",
            po::value<double>(&m_max_D_bias)->default_value(0),
            "Max change in bias per iteration (units)")(

            "equil_steps",
            po::value<long long int>(&m_equil_steps)->default_value(0),
            "Number of equilibration steps")(

            "max_equil_dur",
            po::value<long long int>(&m_max_equil_dur)->default_value(0),
            "Maximum duration of equilibration (s)")(

            "iter_steps",
            po::value<long long int>(&m_iter_steps)->default_value(0),
            "Number of steps per iteration")(

            "iter_swaps",
            po::value<long long int>(&m_iter_swaps)->default_value(0),
            "Maximum number of swaps per iteration")(

            "max_iter_dur",
            po::value<long long int>(&m_max_iter_dur)->default_value(0),
            "Maximum duration of each iteration (s)")(

            "max_rel_P_diff",
            po::value<double>(&m_max_rel_P_diff)->default_value(0.1),
            "Maximum allowed change in P for convergence")(

            "read_biases",
            po::value<bool>(&m_read_biases)->default_value(false),
            "Read US biases from file")(

            "biases_file",
            po::value<string>(&m_biases_file)->default_value(""),
            "Initial guesses at grid biases")(

            "biases_filebase",
            po::value<string>(&m_biases_filebase)->default_value(""),
            "Filebase for grid bias files")(

            "multi_window",
            po::value<bool>(&m_multi_window)->default_value(false),
            "Use multiple windows")(

            "windows_file",
            po::value<string>(&m_windows_file)->default_value(""),
            "File containing windows as min/max pairs of tuples");

    displayed_options.add(us_options);

    po::options_description out_options {"Output options"};
    out_options.add_options()(
            "output_filebase",
            po::value<string>(&m_output_filebase)->default_value(""),
            "Base name for output files")(

            "logging_freq",
            po::value<int>(&m_logging_freq)->default_value(0),
            "Logging frequency")(

            "configs_output_freq",
            po::value<int>(&m_configs_output_freq)->default_value(0),
            "Configuration output write frequency")(

            "vtf_output_freq",
            po::value<int>(&m_vtf_output_freq)->default_value(0),
            "Configuration output write frequency")(

            "vcf_per_domain",
            po::value<bool>(&m_vcf_per_domain)->default_value(false),
            "Write a VCF entry for every domain grown")(

            "counts_output_freq",
            po::value<int>(&m_counts_output_freq)->default_value(0),
            "Counts output write frequency")(

            "times_output_freq",
            po::value<int>(&m_times_output_freq)->default_value(0),
            "Step timing output write frequency")(

            "energies_output_freq",
            po::value<int>(&m_energies_output_freq)->default_value(0),
            "Energies output write frequency")(

            "ops_to_output",
            po::value<string>(),
            "Order parameters to output to file")(

            "order_params_output_freq",
            po::value<int>(&m_order_params_output_freq)->default_value(0),
            "Order parameters write frequency")(

            "rand_engine_state_output_freq",
            po::value<int>(&m_rand_engine_state_output_freq)->default_value(0),
            "Random engine state write frequency")(

            "vmd_pipe_freq",
            po::value<int>(&m_vmd_pipe_freq)->default_value(0),
            "Realtime VMD visualization updating frequency")(

            "create_vmd_instance",
            po::value<bool>(&m_create_vmd_instance)->default_value(false),
            "Create VMD instance");

    displayed_options.add(out_options);

    // Command line variable map
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, displayed_options), vm);
    po::notify(vm);
    if (vm.count("help")) {
        cout << std::endl;
        cout << displayed_options;
        cout << std::endl;
        exit(1);
    }
    if (vm.count("version")) {
        cout << std::endl;
        cout << "Git commit hash: " << GIT_COMMIT << std::endl;
        cout << std::endl;
        exit(1);
    }
    if (not vm.count("parameter_filename")) {
        cout << "Input parameter file must be provided" << std::endl
             << "Run with -h to see all options" << std::endl;
        cout << std::endl;
        exit(1);
    }

    // Parameter file variable map
    ifstream parameter_file {vm["parameter_filename"].as<string>()};
    po::store(po::parse_config_file(parameter_file, displayed_options), vm);
    po::notify(vm);

    process_custom_types(vm);
}

void InputParameters::process_custom_types(po::variables_map vm) {
    if (vm.count("excluded_staples")) {
        string excluded_staples_s {vm["excluded_staples"].as<string>()};
        m_excluded_staples = string_to_int_vector(excluded_staples_s);
    }
    if (vm.count("restart_traj_files")) {
        string file_s = vm["restart_traj_files"].as<string>();
        m_restart_traj_files = string_to_string_vector(file_s);
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
    if (vm.count("stacking_mults")) {
        string stacking_mults_s {vm["stacking_mults"].as<string>()};
        m_stacking_mults = string_to_double_vector(stacking_mults_s);
    }
    if (vm.count("ops_to_output")) {
        string ops_to_output_s {vm["ops_to_output"].as<string>()};
        if (ops_to_output_s != "") {
            m_ops_to_output = string_to_string_vector(ops_to_output_s);
        }
    }
}
} // namespace parser
