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
                    int max_num_staples,
                    OrigamiSystem& origami_system);
            virtual ~OrigamiOutputFile() {};

            void open_write_close();
            virtual void write(long int step) = 0;

            const string m_filename;
            const int m_write_freq;

        protected:
            OrigamiSystem& m_origami_system;
            ofstream m_file;
            int m_max_num_domains;
    };

    class OrigamiVSFOutputFile: public OrigamiOutputFile {
        // Write a VSF file of the system
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int);
    };

    class OrigamiTrajOutputFile: public OrigamiOutputFile {
        // Trajectory output file for simulation configurations
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step);
    };

    class OrigamiVCFOutputFile: public OrigamiOutputFile {
        // Output file for simulation configurations compatible with VMD
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step);
    };

    class OrigamiOrientationOutputFile: public OrigamiOutputFile {
        // Simple format for orientation vector output mainly for vmd viewing
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step);
    };

    class OrigamiStateOutputFile: public OrigamiOutputFile {
        // Outputs binding state of each domain
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step);
    };

    class OrigamiEnergiesOutputFile: public OrigamiOutputFile {
        public:
            OrigamiEnergiesOutputFile(
                    string filename,
                    int write_freq,
                    int max_num_staples,
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
                    int max_num_staples,
                    OrigamiSystem& origami_system);
            void write(long int step);
        private:
            SystemOrderParams* m_system_order_params;
            vector<OrderParam*> m_order_params;
    };

}

#endif // FILES_H
