// files.cpp

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "json/json.h"

#include "files.h"
#include "order_params.h"

using std::ifstream;
using std::vector;
using std::cout;

using namespace Files;
using namespace OrderParams;

OrigamiInputFile::OrigamiInputFile(string filename) {
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

OrigamiOutputFile::OrigamiOutputFile(
        string filename,
        int write_freq,
        OrigamiSystem& origami_system) :
        m_filename {filename},
        m_write_freq {write_freq},
        m_origami_system {origami_system} {

    m_file.open(m_filename);
}

void OrigamiTrajOutputFile::write(long int step) {
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

void OrigamiCountsOutputFile::write(long int step) {
    m_file << step << " ";
    m_file << m_origami_system.num_staples() << " ";
    m_file << m_origami_system.num_unique_staples() << " ";
    m_file << m_origami_system.num_bound_domain_pairs() << " ";
    m_file << m_origami_system.num_fully_bound_domain_pairs() << " ";
    m_file << m_origami_system.num_misbound_domain_pairs() << " ";
    m_file << "\n";
}

OrigamiEnergiesOutputFile::OrigamiEnergiesOutputFile(string filename,
        int write_freq, OrigamiSystem& origami_system) :
        OrigamiOutputFile {filename, write_freq, origami_system} {

    m_file << "step energy bias\n";
}

void OrigamiEnergiesOutputFile::write(long int step) {
    m_file << step << " ";
    m_file << m_origami_system.energy() << " ";
    m_file << m_origami_system.bias() << " ";
    m_file << "\n";
}

OrigamiOrderParamsOutputFile::OrigamiOrderParamsOutputFile(string filename,
        int write_freq, OrigamiSystem& origami_system) :
        OrigamiOutputFile {filename, write_freq, origami_system},
        m_system_order_params {origami_system.get_system_order_params()},
        m_order_params {m_system_order_params->get_order_params()} {
    for (auto order_param: m_order_params) {
        m_file << order_param->get_label() << " ";
    }
    m_file << "\n";
}

void OrigamiOrderParamsOutputFile::write(long int step) {
    m_file << step;
    for (auto order_param: m_order_params) {
        m_file << " " << order_param->calc_param();
    }
    m_file << "\n";
}
