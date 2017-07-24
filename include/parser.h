// parser.h

#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <vector>

namespace Parser {

    using std::string;
    using std::vector;

    enum class MovetypeID {
        OrientationRotation,
        MetStapleExchange,
        MetStapleRegrowth,
        CBStapleRegrowth,
        CTCBScaffoldRegrowth,
        CTCBLinkerRegrowth
    };

    // String parsing helper functions
    vector<int> string_to_int_vector(string string_v);
    vector<double> string_to_double_vector(string string_v);
    vector<string> string_to_string_vector(string string_v);

    // These are from a stackexchange answer
    template<typename Out>
    void split(const string& s, char delim, Out result);

    vector<string> split(const string& s, char delim);

    class Fraction {
        public:
            Fraction(string unparsed_fraction);
            inline double to_double() const {return m_double_fraction;}
        private:
            double m_double_fraction;
            double m_numerator;
            double m_denominator;
    };

    class InputParameters {
        public:
            InputParameters(int argc, char* argv[]);

            // This constructor is for testing purposes
            InputParameters() {};

            // System input parameters
            string m_origami_input_filename;
            double m_temp; // K
            double m_staple_M; // mol/L
            double m_cation_M; // mol/L
            double m_temp_for_staple_u; // K
            double m_staple_u_mult;
            double m_lattice_site_volume; // L
            bool m_cyclic;
            bool m_no_misbinding;
            string m_energy_filebase;
            string m_restart_traj_file;
            vector<string> m_restart_traj_files {};
            string m_restart_traj_filebase;
            string m_restart_traj_postfix;
            int m_restart_step;
            vector<int> m_restart_steps;
            string m_vmd_file_dir;

            // Order parameters
            bool m_distance_sum;
            vector<int> m_distance_pairs {};
            bool m_grid_bias;

            // Bias functions
            bool m_biases_present {false};
            bool m_distance_bias;
            int m_min_dist;
            int m_max_dist;
            double m_max_bias;
            double m_bias_mult;
            bool m_square_well_bias;
            int m_max_well_param;
            int m_min_well_param;
            double m_well_bias;
            double m_outside_bias;

            // General simulation parameters
            string m_simulation_type;
            long long int m_steps;
            int m_logging_freq;
            int m_centering_freq;
            int m_centering_domain;
            int m_constraint_check_freq;
            vector<MovetypeID> m_movetypes {};
            vector<double> m_movetype_probs {};
            int m_max_displacement;
            int m_max_turns;
            string m_num_walks_filename;
            int m_exchange_mult;
            int m_max_total_staples;
            int m_max_type_staples;

            // Enumerator parameters
            bool m_enumerate_staples_only;

            // Annealing simulation parameters
            double m_max_temp;
            double m_min_temp;
            double m_temp_interval;
            int m_steps_per_temp;

            // Parallel tempering simulation parameters
            vector<double> m_temps {};
            int m_num_reps;
            int m_exchange_interval;
            bool m_constant_staple_M;
            vector<double> m_bias_mults {};
            vector<double> m_chem_pot_mults {};
            string m_restart_swap_file;

            // Umbrella sampling simulation parameters
            vector<string> m_order_params {};
            int m_max_num_iters;
            double m_max_D_bias;
            long long int m_equil_steps;
            long long int m_iter_steps;
            long long int m_prod_steps;
            double m_max_rel_P_diff;
            string m_biases_file;
            string m_biases_filebase;
            bool m_multi_window;
            string m_windows_file;

            // Output options
            string m_output_filebase;
            int m_configs_output_freq;
            int m_vtf_output_freq;
            bool m_vcf_per_domain;
            int m_counts_output_freq;
            int m_energies_output_freq;
            int m_order_params_output_freq;
            int m_vmd_pipe_freq;
            bool m_create_vmd_instance;

        private:
            void set_default_sim_options();
    };
}

#endif // PARSER_H
