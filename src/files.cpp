// files.cpp

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "files.h"
#include "utility.h"

namespace files {

    using std::ifstream;
    using std::vector;
    using std::cout;

    using utility::FileMisuse;
    using utility::Occupancy;
    using utility::VectorThree;

    OrigamiInputFile::OrigamiInputFile(string filename) {
        try {
            read_file(filename);
        }
        catch (Json::RuntimeError) {
            cout << "Problem reading origami system file\n";
            throw;
        }
    }

    void OrigamiInputFile::read_file(string filename) {
        ifstream jsonraw {filename, ifstream::binary};
        Json::Value jsonroot;
        jsonraw >> jsonroot;

        // Extract sequences
        Json::Value jsonseqs {jsonroot["origami"]["sequences"]};
        for (unsigned int i {0}; i != jsonseqs.size(); i++) {
            m_sequences.push_back({});
            for (unsigned int j {0}; j != jsonseqs[i].size(); j++) {
                m_sequences[i].push_back(jsonseqs[i][j].asString());
            }
        }

        // Extract identities
        Json::Value jsonidents {jsonroot["origami"]["identities"]};
        for (unsigned int i {0}; i != jsonidents.size(); i++) {
            m_identities.push_back({});
            auto num_domains {jsonidents[i].size()};
            for (unsigned int j {0}; j != num_domains; j++) {
                int ident {jsonidents[i][j].asInt()};
                m_identities[i].push_back(ident);
            }
        }

        // Extract configuration
        // Note for now always using first configuration, may change file format later
        Json::Value jsonconfig {jsonroot["origami"]["configurations"][0]["chains"]};
        for (unsigned int i {0}; i != jsonconfig.size(); i++) {
            Json::Value jsonchain {jsonconfig[i]};
            int index {jsonchain["index"].asInt()};
            int identity {jsonchain["identity"].asInt()};
            vector<VectorThree> positions {};
            vector<VectorThree> orientations {};
            for (unsigned int j {0}; j != jsonchain["positions"].size(); j++) {
                int posx {jsonchain["positions"][j][0].asInt()};
                int posy {jsonchain["positions"][j][1].asInt()};
                int posz {jsonchain["positions"][j][2].asInt()};
                positions.push_back(VectorThree {posx, posy, posz});
                int orex {jsonchain["orientations"][j][0].asInt()};
                int orey {jsonchain["orientations"][j][1].asInt()};
                int orez {jsonchain["orientations"][j][2].asInt()};
                orientations.push_back(VectorThree {orex, orey, orez});
            }
            Chain chain {index, identity, positions, orientations};
            m_chains.push_back(chain);
        }

        // Extract cyclic flag
        Json::Value jsoncyclic {jsonroot["origami"]["cyclic"]};
        m_cyclic = jsoncyclic.asBool();
    }

    vector<vector<int>> OrigamiInputFile::get_identities() {
        return m_identities;
    }

    vector<vector<string>> OrigamiInputFile::get_sequences() {
        return m_sequences;
    }

    vector<Chain> OrigamiInputFile::get_config() {
        return m_chains;
    }

    bool OrigamiInputFile::is_cyclic() {
        return m_cyclic;
    }

    OrigamiTrajInputFile::OrigamiTrajInputFile(string filename):
            m_filename {filename} {
        m_file.open(m_filename);
    }
    
    vector<Chain> OrigamiTrajInputFile::read_config(int step) {
        try {
            return internal_read_config(step);
        }
        catch (Json::RuntimeError) {
            cout << "Problem reading trajectory file config\n";
            throw;
        }
    }

    vector<Chain> OrigamiTrajInputFile::internal_read_config(int step) {
        vector<Chain> step_chains {};
        go_to_step(step);
        while (true) {
            string identity_line;
            std::getline(m_file, identity_line);
            if (identity_line.empty()) {
                break;
            }

            std::istringstream identity_line_stream {identity_line};
            int chain_index;
            identity_line_stream >> chain_index;
            int chain_identity;
            identity_line_stream >> chain_identity;

            string pos_line;
            std::getline(m_file, pos_line);
            std::istringstream pos_line_stream {pos_line};
            vector<VectorThree> positions {};
            while (not pos_line_stream.eof()) {
                int x;
                pos_line_stream >> x;
                int y;
                pos_line_stream >> y;
                int z;
                pos_line_stream >> z;
                positions.push_back({x, y, z});
            }

            string ore_line;
            std::getline(m_file, pos_line);
            std::istringstream ore_line_stream {pos_line};
            vector<VectorThree> orientations {};
            while (not ore_line_stream.eof()) {
                int x;
                ore_line_stream >> x;
                int y;
                ore_line_stream >> y;
                int z;
                ore_line_stream >> z;
                orientations.push_back({x, y, z});
            }

            Chain chain {chain_index, chain_identity, positions, orientations};
            step_chains.push_back(chain);

        }

        return step_chains;
    }

    void OrigamiTrajInputFile::go_to_step(unsigned int step){
        // Returns line after step number
        // Really ugly fragile method for doing this
        m_file.seekg(std::ios::beg);
        for(unsigned int i=0; i != step; ++i){
            bool end_of_step_reached {false};
            while (not end_of_step_reached) {
                string line;
                std::getline(m_file, line);
                if (line.empty()) {
                    end_of_step_reached = true;
                }
                if (m_file.eof()) {
                    throw FileMisuse {};
                }
            }
        }

        // Check that read step is requested step
        string step_s;
        std::getline(m_file, step_s);
        //int read_step {std::stoi(step_s)};
        //if (read_step != step) {
        //    throw FileMisuse {};
        //}
    }

    OrigamiMovetypeFile::OrigamiMovetypeFile(string filename):
            m_filename {filename} {

        try {
            read_file(filename);
        }
        catch (Json::RuntimeError) {
            cout << "Problem reading movetype file\n";
            throw;
        }
    }

    void OrigamiMovetypeFile::read_file(string filename) {
        ifstream jsonraw {filename, ifstream::binary};
        Json::Value jsonroot;
        jsonraw >> jsonroot;
        m_jsonmovetypes = jsonroot["origami"]["movetypes"];
        for (unsigned int i {0}; i != m_jsonmovetypes.size(); i++) {
            Json::Value jsonmovetype {m_jsonmovetypes[i]};
            string type {jsonmovetype["type"].asString()};
            string label {jsonmovetype["label"].asString()};
            string freq_raw {jsonmovetype["freq"].asString()};
            utility::Fraction freq_frac {freq_raw};
            double freq {freq_frac.to_double()};
            m_types.push_back(type);
            m_labels.push_back(label);
            m_freqs.push_back(freq);
        }
    }

    vector<string> OrigamiMovetypeFile::get_types() {
        return m_types;
    }

    vector<string> OrigamiMovetypeFile::get_labels() {
        return m_labels;
    }

    vector<double> OrigamiMovetypeFile::get_freqs() {
        return m_freqs;
    }

    bool OrigamiMovetypeFile::get_bool_option(int movetype_i, string key) {
        return m_jsonmovetypes[movetype_i][key].asBool();
    }
     
    double OrigamiMovetypeFile::get_double_option(int movetype_i, string key) {
        return m_jsonmovetypes[movetype_i][key].asDouble();
    }
     
    vector<double> OrigamiMovetypeFile::get_double_vector_option(int movetype_i,
            string key) {
        vector<double> dv {};
        for (unsigned int i {0}; i != m_jsonmovetypes[movetype_i][key].size(); i++) {
            dv.push_back(m_jsonmovetypes[movetype_i][key][i].asDouble());
        }

        return dv;
    }

    string OrigamiMovetypeFile::get_string_option(int movetype_i, string key) {
        return m_jsonmovetypes[movetype_i][key].asString();
    }

    int OrigamiMovetypeFile::get_int_option(int movetype_i, string key) {
        return m_jsonmovetypes[movetype_i][key].asInt();
    }

    vector<vector<string>> OrigamiLeveledInput::get_types_by_level() {
        return m_level_to_types;
    }

    vector<vector<string>> OrigamiLeveledInput::get_labels_by_level() {
        return m_level_to_labels;
    }

    vector<vector<string>> OrigamiLeveledInput::get_tags_by_level() {
        return m_level_to_tags;
    }

    double OrigamiLeveledInput::get_double_option(int i, int j, string key) {
        int json_i {m_level_to_indices[i][j]};
        return m_json_ops[json_i][key].asDouble();
    }

    string OrigamiLeveledInput::get_string_option(int i, int j, string key) {
        int json_i {m_level_to_indices[i][j]};
        return m_json_ops[json_i][key].asString();
    }

    int OrigamiLeveledInput::get_int_option(int i, int j, string key) {
        int json_i {m_level_to_indices[i][j]};
        return m_json_ops[json_i][key].asInt();
    }

    bool OrigamiLeveledInput::get_bool_option(int i, int j, string key) {
        int json_i {m_level_to_indices[i][j]};
        return m_json_ops[json_i][key].asBool();
    }

    vector<string> OrigamiLeveledInput::get_vector_string_option(int i, int j, string key) {
        int json_i {m_level_to_indices[i][j]};
        vector<string> sv {};
        for (unsigned int i {0}; i != m_json_ops[json_i][key].size(); i++) {
            sv.push_back(m_json_ops[json_i][key][i].asString());
        }
        return sv;
    }

    OrigamiOrderParamsFile::OrigamiOrderParamsFile(string filename) {
        m_filename  = filename;
        try {
            read_file(filename);
        }
        catch (Json::RuntimeError) {
            cout << "Problem reading order parameter file\n";
            throw;
        }
    }

    void OrigamiOrderParamsFile::read_file(string filename) {
        ifstream jsonraw {filename, ifstream::binary};
        Json::Value jsonroot;
        jsonraw >> jsonroot;
        m_json_ops = jsonroot["origami"]["order_params"];
        int max_level {0};
        vector<string> types {};
        vector<string> labels {};
        vector<string> tags {};
        vector<int> levels {};
        for (unsigned int i {0}; i != m_json_ops.size(); i++) {
            Json::Value json_op {m_json_ops[i]};
            string type {json_op["type"].asString()};
            string label {json_op["label"].asString()};
            string tag {json_op["tag"].asString()};
            int level {json_op["level"].asInt()};
            if (level > max_level) {
                max_level = level;
            }
            types.push_back(type);
            labels.push_back(label);
            tags.push_back(tag);
            levels.push_back(level);
        }

        for (int i {0}; i != max_level + 1; i++) {
            m_level_to_types.push_back({});
            m_level_to_labels.push_back({});
            m_level_to_tags.push_back({});
            m_level_to_indices.push_back({});
        }
        for (unsigned int i {0}; i != m_json_ops.size(); i++) {
            m_level_to_types[levels[i]].push_back(types[i]);
            m_level_to_labels[levels[i]].push_back(labels[i]);
            m_level_to_tags[levels[i]].push_back(tags[i]);
            m_level_to_indices[levels[i]].push_back(i);
        }
    }

    OrigamiBiasFunctionsFile::OrigamiBiasFunctionsFile(string filename) {
        m_filename  = filename;
        try {
            read_file(filename);
        }
        catch (Json::RuntimeError) {
            cout << "Problem reading bias function file\n";
            throw;
        }
    }

    void OrigamiBiasFunctionsFile::read_file(string filename) {
        ifstream jsonraw {filename, ifstream::binary};
        Json::Value jsonroot;
        jsonraw >> jsonroot;
        m_json_ops = jsonroot["origami"]["bias_functions"];
        int max_level {0};
        vector<string> types {};
        vector<string> labels {};
        vector<string> tags {};
        vector<int> levels {};
        vector<vector<string>> op_tags {};
        vector<vector<string>> d_biases_tags {};
        for (unsigned int i {0}; i != m_json_ops.size(); i++) {
            Json::Value json_op {m_json_ops[i]};
            string type {json_op["type"].asString()};
            string label {json_op["label"].asString()};
            string tag {json_op["tag"].asString()};
            int level {json_op["level"].asInt()};
            vector<string> ops {};
            for (unsigned int j {0}; j != m_json_ops[i]["ops"].size(); j++) {
                ops.push_back(m_json_ops[i]["ops"][j].asString());
            }
            string d_biases_raw {m_json_ops[i]["bias_funcs"].asString()};
            if (level > max_level) {
                max_level = level;
            }
            types.push_back(type);
            labels.push_back(label);
            tags.push_back(tag);
            levels.push_back(level);
            op_tags.push_back(ops);
            d_biases_tags.push_back(utility::string_to_string_vector(d_biases_raw));
        }

        for (int i {0}; i != max_level + 1; i++) {
            m_level_to_types.push_back({});
            m_level_to_labels.push_back({});
            m_level_to_tags.push_back({});
            m_level_to_indices.push_back({});
            m_level_to_ops.push_back({});
            m_level_to_d_biases.push_back({});
        }
        for (unsigned int i {0}; i != m_json_ops.size(); i++) {
            m_level_to_types[levels[i]].push_back(types[i]);
            m_level_to_labels[levels[i]].push_back(labels[i]);
            m_level_to_tags[levels[i]].push_back(tags[i]);
            m_level_to_indices[levels[i]].push_back(i);
            m_level_to_ops[levels[i]].push_back(op_tags[i]);
            m_level_to_d_biases[levels[i]].push_back(d_biases_tags[i]);
        }
    }

    vector<vector<vector<string>>> OrigamiBiasFunctionsFile::get_ops_by_level() {
        return m_level_to_ops;
    }

    vector<vector<vector<string>>> OrigamiBiasFunctionsFile::get_d_biases_by_level() {
        return m_level_to_d_biases;
    }

    OrigamiOutputFile::OrigamiOutputFile(
            string filename,
            int write_freq,
            int max_num_staples,
            OrigamiSystem& origami_system) :
            m_filename {filename},
            m_write_freq {write_freq},
            m_origami_system {origami_system} {

        // This assumes 2 domain staples
        m_max_num_domains = 2 * max_num_staples +
                m_origami_system.get_chain(0).size();
        m_file.open(m_filename);
    }

    void OrigamiOutputFile::open_write_close() {
        m_file.open(m_filename);
        write(0, 0);
        m_file.close();
    }

    void OrigamiVSFOutputFile::write(long int, double) {
        m_file << "atom 0:";
        size_t scaffold_length {m_origami_system.get_chain(0).size()};
        m_file << scaffold_length << " radius 0.25 type scaffold\n";

        m_file << "atom " << scaffold_length << ":" << m_max_num_domains - 1;
        m_file << " radius 0.25 type staple";
    }

    void OrigamiTrajOutputFile::write(long int step, double) {
        m_file << step << "\n";
        for (auto chain: m_origami_system.get_chains()) {
            m_file << chain[0]->m_c << " " << chain[0]->m_c_ident << "\n";
            for (auto domain: chain) {
                for (int i {0}; i != 3; i++) {
                    m_file << domain->m_pos[i] << " ";
                }
            }
            m_file << "\n";
            for (auto domain: chain) {
                for (int i {0}; i != 3; i++) {
                    m_file << domain->m_ore[i] << " ";
                }
            }
            m_file << "\n";
        }
        m_file << "\n";
    }

    void OrigamiVCFOutputFile::write(long int, double) {
        m_file << "timestep\n";
        for (auto chain: m_origami_system.get_chains()) {
            for (auto domain: chain) {
                for (int i {0}; i != 3; i++) {
                    m_file << domain->m_pos[i] << " ";
                }
                m_file << "\n";
            }
        }

        int num_doms {m_origami_system.num_domains()};
        for (int dom_i {num_doms}; dom_i != m_max_num_domains; dom_i++) {
            for (int i {0}; i !=3; i++) {
                m_file << 0 << " ";
            }
            m_file << "\n";
        }
        m_file << "\n";
    }

    void OrigamiOrientationOutputFile::write(long int, double) {
        for (auto chain: m_origami_system.get_chains()) {
            for (auto domain: chain) {
                for (int i {0}; i != 3; i++) {
                    if (domain->m_state != Occupancy::unassigned) {
                        m_file << domain->m_ore[i] << " ";
                    }
                    else {
                        m_file << 0 << " ";
                    }
                }
            }
        }
        int num_doms {m_origami_system.num_domains()};
        for (int dom_i {num_doms}; dom_i != m_max_num_domains; dom_i++) {
            for (int i {0}; i !=3; i++) {
                m_file << 0 << " ";
            }
        }
        m_file << "\n";
    }

    void OrigamiStateOutputFile::write(long int, double) {
        for (auto chain: m_origami_system.get_chains()) {
            for (auto domain: chain) {
                if (domain->m_state == Occupancy::unbound) {
                    m_file << "1 ";
                }
                else if (domain->m_state == Occupancy::bound) {
                    m_file << "2 ";
                }
                else if (domain->m_state == Occupancy::misbound) {
                    m_file << "3 ";
                }
                else {
                    m_file << "0 ";
                }
            }
        }
        int num_doms {m_origami_system.num_domains()};
        for (int dom_i {num_doms}; dom_i != m_max_num_domains; dom_i++) {
            m_file << "-1 ";
        }
                
        m_file << "\n";
    }

    void OrigamiCountsOutputFile::write(long int step, double) {
        m_file << step << " ";
        m_file << m_origami_system.num_staples() << " ";
        m_file << m_origami_system.num_unique_staples() << " ";
        m_file << m_origami_system.num_bound_domain_pairs() << " ";
        m_file << m_origami_system.num_fully_bound_domain_pairs() << " ";
        m_file << m_origami_system.num_misbound_domain_pairs() << " ";
        m_file << "\n";
    }

    void OrigamiStaplesBoundOutputFile::write(long int step, double) {
        m_file << step << " ";
        for (auto staple_count: m_origami_system.get_staple_counts()) {
            m_file << staple_count << " ";
        }
        m_file << "\n";
    }

    void OrigamiStaplesFullyBoundOutputFile::write(long int step, double) {
        m_file << step << " ";
        for (size_t staple_ident {1}; staple_ident !=
                m_origami_system.m_identities.size(); staple_ident++) {
            int staple_state {0};
            for (auto staple_i: m_origami_system.staples_of_ident(staple_ident)) {
                bool fully_bound {true};
                for (auto domain: m_origami_system.get_chain(staple_i)) {
                    if (domain->m_state != Occupancy::bound) {
                        fully_bound = false;
                        break;
                    }
                }
                if (fully_bound) {
                    staple_state = 1;
                    break;
                }
            }
            m_file << staple_state << " ";
        }
        m_file << "\n";
    }

    OrigamiEnergiesOutputFile::OrigamiEnergiesOutputFile(
            string filename,
            int write_freq,
            int max_num_staples,
            OrigamiSystem& origami_system,
            SystemBiases& biases) :
            OrigamiOutputFile {filename, write_freq, max_num_staples, origami_system},
            m_biases {biases} {

        m_file << "step tenergy henthalpy hentropy stacking bias\n";
        m_file.precision(10);
    }

    void OrigamiEnergiesOutputFile::write(long int step, double) {
        m_file << step << " ";
        m_file << m_origami_system.energy() << " ";
        m_origami_system.update_enthalpy_and_entropy();
        m_file << m_origami_system.hybridization_enthalpy() << " ";
        m_file << m_origami_system.hybridization_entropy() << " ";
        m_file << m_origami_system.stacking_energy() << " ";
        m_file << m_biases.get_total_bias() << " ";
        m_file << "\n";
    }

    OrigamiTimesOutputFile::OrigamiTimesOutputFile(
            string filename,
            int write_freq,
            int max_num_staples,
            OrigamiSystem& origami_system):
            OrigamiOutputFile {filename, write_freq, max_num_staples, origami_system} {

        m_file << "step time\n";
    }

    void OrigamiTimesOutputFile::write(long int step, double time) {
        m_file << step << " ";
        m_file << time;
        m_file << "\n";
    }

    OrigamiOrderParamsOutputFile::OrigamiOrderParamsOutputFile(
            string filename,
            int write_freq,
            int max_num_staples,
            OrigamiSystem& origami_system,
            SystemOrderParams& ops,
            vector<string> op_tags) :
            OrigamiOutputFile {filename, write_freq, max_num_staples, origami_system},
            m_ops {ops} {

        for (auto op_tag: op_tags) {
            OrderParam& op {m_ops.get_order_param(op_tag)};
            m_file << op_tag << ", ";
            m_file.precision(10);
            m_ops_to_output.emplace_back(op);
        }
        m_file << "\n";
    }

    void OrigamiOrderParamsOutputFile::write(long int step, double) {
        m_file << step;
        for (auto op: m_ops_to_output) {
            m_file << " " << op.get().calc_param();
        }
        m_file << "\n";
    }
}
