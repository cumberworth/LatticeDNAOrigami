// parser.h

#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace Parser {
    vector<double> string_to_double_vector(string string_v);

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

            // System input parameters
            string m_origami_input_filename;

            // Kelvin
            double m_temp {300};

            // mol/L
            double m_staple_M {1};

            // mol/L
            double m_cation_M {1};

            // Kelvin
            double m_temp_for_staple_u {m_temp};

            // L
            double m_lattice_site_volume {1};
            bool m_cyclic {false};
            string m_energy_filebase {""};
            string m_restart_traj_file {""};
            int m_restart_step;

            // General simulation parameters
            string m_simulation_type {"constant_temp"};
            long int m_steps {0};
            int m_logging_freq {0};
            int m_centering_freq {0};
            vector<int> m_movetypes {};
            vector<double> m_movetype_probs {};
            string m_num_walks_filename {""};
            int m_exchange_mult {1};

            // Annealing simulation parameters
            double m_max_temp {};
            double m_min_temp {};
            double m_temp_interval {};
            int m_steps_per_temp {};

            // Parallel tempring simulation parameters
            vector<double> m_temps {};
            int m_num_reps {};
            int m_exchange_interval {};
            bool m_constant_staple_M {true};

            // Output options
            string m_output_filebase;
            int m_configs_output_freq {0};
            int m_counts_output_freq {0};

    };
}

#endif // PARSER_H
