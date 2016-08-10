// files.h

#ifndef FILES_H
#define FILES_H

#include <string>
#include <vector>

#include "origami_system.h"

using namespace Origami;

namespace Files {

    class OrigamiSystemInputFile {
        // Input file for OrigamiSystem configuration and topology
        public:
            OrigamiSystemInputFile(string filename);

            // Properties
            vector<vector<int>> m_identities;
            vector<vector<string>> m_sequences;
            vector<Chain> m_chains;
    };

    class OrigamiSystemOutputFile {
        // Output file for OrigamiSystem configuration and topology
    };

    class OrigamiSystemTrajOutputFile {
        // Trajectory output file for simulation configurations
    };
}

#endif // FILES_H
