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
            OrigamiOutputFile(
                    string filename,
                    int write_freq,
                    OrigamiSystem& origami_system);

            virtual void write(int step) = 0;

            const string m_filename;
            const int m_write_freq;

        protected:
            OrigamiSystem& m_origami_system;
            ofstream m_file;
    };

    class OrigamiTrajOutputFile: public OrigamiOutputFile {
        // Trajectory output file for simulation configurations
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(int step);
    };

    class OrigamiCountsOutputFile: public OrigamiOutputFile {
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(int step);
    };

}

#endif // FILES_H
