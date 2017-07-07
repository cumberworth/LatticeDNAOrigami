// files.h

#ifndef FILES_H
#define FILES_H

#include <string>
#include <vector>
#include <fstream>

#include "origami_system.h"
#include "order_params.h"

namespace Files {

    using std::ofstream;
    using std::ifstream;

    using namespace Origami;
    using namespace OrderParams;

    class OrigamiInputFile {
        // Input file for OrigamiSystem configuration and topology
        public:
            OrigamiInputFile(string filename);
            ~OrigamiInputFile() {};

            vector<vector<int>> get_identities();
            vector<vector<string>> get_sequences();
            vector<Chain> get_config();
            bool is_cyclic();

        private:
            vector<vector<int>> m_identities;
            vector<vector<string>> m_sequences;
            vector<Chain> m_chains;
            bool m_cyclic;
    };

    class OrigamiTrajInputFile {
        // Reading configurations from trajectories
        public:
            OrigamiTrajInputFile(string filename);

            vector<Chain> read_config(int step);

        private:
            string m_filename;
            ifstream m_file;

            void go_to_step(unsigned int num);

    };

    class OrigamiOutputFile {
        // Output file interface
        public:
            OrigamiOutputFile(
                    string filename,
                    int write_freq,
                    OrigamiSystem& origami_system);
            virtual ~OrigamiOutputFile() {};

            virtual void write(long int step) = 0;

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
            void write(long int step);
    };

    class OrigamiEnergiesOutputFile: public OrigamiOutputFile {
        public:
            OrigamiEnergiesOutputFile(
                    string filename,
                    int write_freq,
                    OrigamiSystem& origami_system);
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step);
    };

    class OrigamiCountsOutputFile: public OrigamiOutputFile {
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step);
    };

    class OrigamiOrderParamsOutputFile: public  OrigamiOutputFile {
        public:
            OrigamiOrderParamsOutputFile(
                    string filename,
                    int write_freq,
                    OrigamiSystem& origami_system);
            void write(long int step);
        private:
            SystemOrderParams* m_system_order_params;
            vector<OrderParam*> m_order_params;
    };

}

#endif // FILES_H
