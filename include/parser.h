// parse.h

#ifndef PARSER_H
#define PARSER_H

#include <string>

#include "movetypes.h"

using std::string;

using namespace Movetypes;

namespace Parser {
    class InputParameters {
        public:
            InputParameters(int argc, char* argv[]);

            // Parameters
            string m_origami_input_filename;
            string m_configs_output_filename;
            int m_configs_output_freq {0};

            // Kelvin
            double m_temp {300};

            // mol/L
            double m_staple_M {1};

            // mol/L
            double m_cation_M;

            // L
            double m_lattice_site_volume {1};
            bool m_cyclic {false};
            int m_steps {0};
            int m_logging_freq {0};
            int m_centering_freq {0};
            vector<double> m_movetype_probs {};
            vector<MovetypeConstructor> m_movetype_constructors {};
    };
}

#endif // PARSER_H
