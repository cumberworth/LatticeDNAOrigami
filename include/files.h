// files.h

#ifndef FILES_H
#define FILES_H

#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "json/json.h"

#include "bias_functions.h"
#include "origami_system.h"
#include "order_params.h"

namespace files {

    using std::ifstream;
    using std::ofstream;
    using std::pair;
    using std::reference_wrapper;
    using std::string;
    using std::vector;

    using biasFunctions::SystemBiases;
    using origami::Chain;
    using origami::OrigamiSystem;
    using orderParams::OrderParam;
    using orderParams::SystemOrderParams;

    class OrigamiInputFile {
        // Input file for OrigamiSystem configuration and topology
        public:
            OrigamiInputFile(string filename);
            ~OrigamiInputFile() {};

            vector<vector<int>> get_identities();
            vector<vector<string>> get_sequences();
            vector<double> get_enthalpies();
            vector<double> get_entropies();
            vector<Chain> get_config();
            bool is_cyclic();

        private:
            void read_file(string filename);

            vector<vector<int>> m_identities;
            vector<vector<string>> m_sequences;
            vector<double> m_enthalpies;
            vector<double> m_entropies;
            vector<Chain> m_chains;
            bool m_cyclic;
    };

    class OrigamiTrajInputFile {
        // Reading configurations from trajectories
        public:
            OrigamiTrajInputFile(string filename);

            vector<Chain> read_config(int step);

        private:
            vector<Chain> internal_read_config(int step);
            string m_filename;
            ifstream m_file;

            void go_to_step(unsigned int num);

    };

    class OrigamiMovetypeFile {
        public:
            OrigamiMovetypeFile(string filename);
            vector<string> get_types();
            vector<string> get_labels();
            vector<double> get_freqs();
            bool get_bool_option(int movetype_i, string key);
            double get_double_option(int movetype_i, string key);
            vector<double> get_double_vector_option(int movetype_i, string key);
            string get_string_option(int movetype_i, string key);
            int get_int_option(int movetype_i, string key);

        private:
            void read_file(string filename);

            string m_filename;
            Json::Value m_jsonmovetypes {};

            vector<string> m_types {};
            vector<string> m_labels {};
            vector<double> m_freqs {};
    };

    class OrigamiLeveledInput {
        public:
            virtual ~OrigamiLeveledInput() {}

            vector<vector<string>> get_types_by_level();
            vector<vector<string>> get_tags_by_level();
            vector<vector<string>> get_labels_by_level();

            double get_double_option(int i, int j, string key);
            string get_string_option(int i, int j, string key);
            int get_int_option(int i, int j, string key);
            bool get_bool_option(int i, int j, string key);
            vector<string> get_vector_string_option(int i, int j, string key);

        protected:
            string m_filename;
            Json::Value m_json_ops {};

            vector<vector<string>> m_level_to_types {};
            vector<vector<string>> m_level_to_labels {};
            vector<vector<string>> m_level_to_tags {};
            vector<vector<int>> m_level_to_indices {};
    };

    class OrigamiOrderParamsFile: public OrigamiLeveledInput {
        public:
            OrigamiOrderParamsFile(string filename);

        private:
            void read_file(string filename);

    };

    class OrigamiBiasFunctionsFile: public OrigamiLeveledInput {
        public:
            // This repeats most of the code in the above class
            // Should try and take a more general approach to these json files
            OrigamiBiasFunctionsFile(string filename);

            vector<vector<vector<string>>> get_ops_by_level();
            vector<vector<vector<string>>> get_d_biases_by_level();

        private:
            void read_file(string filename);

            vector<vector<vector<string>>> m_level_to_ops {};
            vector<vector<vector<string>>> m_level_to_d_biases {};
    };

    class OrigamiOutputFile {
        public:
            OrigamiOutputFile(
                    string filename,
                    int write_freq,
                    int max_num_staples,
                    int max_staple_size,
                    OrigamiSystem& origami_system);
            virtual ~OrigamiOutputFile() {};

            void open_write_close();
            virtual void write(long int step, double) = 0;

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
            void write(long int, double);
    };

    class OrigamiTrajOutputFile: public OrigamiOutputFile {
        // Trajectory output file for simulation configurations
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step, double);
    };

    class OrigamiVCFOutputFile: public OrigamiOutputFile {
        // Output file for simulation configurations compatible with VMD
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step, double);
    };

    class OrigamiOrientationOutputFile: public OrigamiOutputFile {
        // Simple format for orientation vector output mainly for vmd viewing
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step, double);
    };

    class OrigamiStateOutputFile: public OrigamiOutputFile {
        // Outputs binding state of each domain
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step, double);
    };

    class OrigamiEnergiesOutputFile: public OrigamiOutputFile {
        public:
            OrigamiEnergiesOutputFile(
                    string filename,
                    int write_freq,
                    int max_num_staples,
                    int max_staple_size,
                    OrigamiSystem& origami_system,
                    SystemBiases& m_biases);
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step, double);

        private:
            SystemBiases& m_biases;
    };

    class OrigamiTimesOutputFile: public OrigamiOutputFile {
        public:
            OrigamiTimesOutputFile(
                    string filename,
                    int write_freq,
                    int max_num_staples,
                    int max_staple_size,
                    OrigamiSystem& origami_system);
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step, double time);
    };

    class OrigamiCountsOutputFile: public OrigamiOutputFile {
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step, double);
    };

    class OrigamiStaplesBoundOutputFile: public OrigamiOutputFile {
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step, double);
    };

    class OrigamiStaplesFullyBoundOutputFile: public OrigamiOutputFile {
        public:
            using OrigamiOutputFile::OrigamiOutputFile;
            void write(long int step, double);
    };

    class OrigamiOrderParamsOutputFile: public  OrigamiOutputFile {
        public:
            OrigamiOrderParamsOutputFile(
                    string filename,
                    int write_freq,
                    int max_num_staples,
                    int max_staple_size,
                    OrigamiSystem& origami_system,
                    SystemOrderParams& ops,
                    vector<string> op_tags);
            void write(long int step, double);

        private:
            SystemOrderParams& m_ops;
            vector<reference_wrapper<OrderParam>> m_ops_to_output {};
    };

}

#endif // FILES_H
