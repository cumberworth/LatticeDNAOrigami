// utility.h

#ifndef UTILITY_H
#define UTILITY_H

#include <array>
#include <list>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace utility {
    
    using std::array;
    using std::bind;
    using std::pair;
    using std::string;
    using std::unique_ptr;
    using std::unordered_map;
    using std::vector;

    // Find index for given element
    int index(vector<int> container, int element);

    // Exception types
    struct NoElement {};
    struct OrigamiMisuse {};
    struct FileMisuse {};
    struct NotImplemented {};
    struct SimulationMisuse {};

    // Constrants
    const double NA {6.022140857e23};

    enum class Occupancy {
          unassigned,
          unbound,
          bound,
          misbound
      };

    // Vector for Z3 space
    class VectorThree {
        public:
            VectorThree(int x, int y, int z): m_container {{x, y, z}} {};
            VectorThree(): m_container {{0, 0, 0}} {};

            VectorThree operator-();

            VectorThree operator+(const VectorThree& v_2) const;
            VectorThree operator-(const VectorThree& v_2) const;
            bool operator!=(const VectorThree& v_2) const;

            int& operator[](const size_t& i) {return m_container[i];};
            const int& at(const size_t& i) const {return m_container.at(i);};

            VectorThree rotate_half(VectorThree axis);
            VectorThree rotate(VectorThree origin, VectorThree axis, int turns);
            int sum();
            int abssum();
            VectorThree absolute();
            VectorThree sort();

        private:
            array<int, 3> m_container;
		    friend class boost::serialization::access;
            template<typename Archive>
            void serialize(Archive& arch, const unsigned int) {
                arch& m_container;
        	}
    };

    // Hash needs this, ambiguous to also have method
    bool operator==(const VectorThree& v1, const VectorThree& v2);

    struct StapleExchangeTracking {
        bool staple_insertion;
        bool no_staples;
        int staple_type;
    };

    bool operator==(
            const StapleExchangeTracking& t1,
            const StapleExchangeTracking& t2);

    struct StapleRegrowthTracking {
        bool no_staples;
        int staple_type;
    };

    bool operator==(
            const StapleRegrowthTracking& t1,
            const StapleRegrowthTracking& t2);

    struct CTCBScaffoldRegrowthTracking {
        int num_scaffold_domains;
        int num_staples;
    };

    bool operator==(
            const CTCBScaffoldRegrowthTracking& t1,
            const CTCBScaffoldRegrowthTracking& t2);

    struct ScaffoldRGRegrowthTracking {
        int num_scaffold_domains;
    };

    bool operator==(
            const ScaffoldRGRegrowthTracking& t1,
            const ScaffoldRGRegrowthTracking& t2);

    struct CTCBLinkerRegrowthTracking {
        bool central_domains_connected;
        int num_linker_domains;
        int num_linker_staples;
        int num_central_domains;
        int num_central_staples;
        int disp_sum;
        int rot_turns;
    };

    bool operator==(
            const CTCBLinkerRegrowthTracking& t1,
            const CTCBLinkerRegrowthTracking& t2);

    // Base unit vectors
    const VectorThree xhat {1, 0, 0};
    const VectorThree yhat {0, 1, 0};
    const VectorThree zhat {0, 0, 1};

    // All possible unit vectors
    const vector<VectorThree> vectors {
            {1, 0, 0},
            {-1, 0, 0},
            {0, 1, 0},
            {0, -1, 0},
            {0, 0, 1},
            {0, 0, -1}};

    const vector<VectorThree> basis_vectors {xhat, yhat, zhat};

    /** Create all possible pairs of a given set (as a vector) */
    template<typename aT>
    vector<pair<aT, aT>> all_pairs(vector<aT> v) {
        vector<pair<aT, aT>> c {};
        for (auto a: v) {
            for (auto b: v) {
                c.emplace_back(a, b);
            }
        }

        return c;
    }

    template<typename aT>
    /**
      * Create map from pair to index
      *
      * For use alongsie the all_pairs method
      */
    unordered_map<pair<aT, aT>, int> pair_to_index(vector<pair<aT, aT>> v) {
        unordered_map<pair<aT, aT>, int> c {};
        int i {0};
        for (auto a: v) {
            c[a] = i;
            i++;
        }

        return c;
    }

    // String parsing helper functions
    vector<int> string_to_int_vector(string string_v);
    vector<double> string_to_double_vector(string string_v);
    vector<string> string_to_string_vector(string string_v);

    // These are from a stackexchange answer
    template<typename Out>
    void split(const string& s, char delim, Out result);

    vector<string> split(const string& s, char delim);

    class Fraction {
        public:
            Fraction(string unparsed_fraction);

            /// test
            inline double to_double() const {return m_double_fraction;}
        private:
            double m_double_fraction;
            double m_numerator;
            double m_denominator;
    };
}

#endif // UTILITY_H
