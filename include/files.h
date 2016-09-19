// files.h

#ifndef FILES_H
#define FILES_H

#include <string>
#include <vector>
#include <fstream>

#include "origami_system.h"

using std::ofstream;

using namespace Origami;

namespace Files {

    class OrigamiInputFile {
        // Input file for OrigamiSystem configuration and topology
        public:
            OrigamiInputFile(string filename);

            // Properties
            vector<vector<int>> m_identities;
            vector<vector<string>> m_sequences;
            vector<Chain> m_chains;
    };

    class OrigamiOutputFile {
        // Output file interface
        public:
            virtual void write(int step) = 0;
            int m_write_freq;
    };

    class OrigamiTrajOutputFile: public OrigamiOutputFile {
        // Trajectory output file for simulation configurations
        public:
            OrigamiTrajOutputFile(string filename, int write_freq, OrigamiSystem& origami_system);
            void write(int step);

        private:
            string m_filename;
            OrigamiSystem& m_origami_system;
            ofstream m_file;
    };
}

#endif // FILES_H
