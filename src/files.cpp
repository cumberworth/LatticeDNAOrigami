// files.cpp

#include <vector>
#include <fstream>
#include <iostream>
#include <json/json.h>

#include "files.h"

using std::ifstream;
using std::vector;
using std::cout;

using namespace Files;

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
}
