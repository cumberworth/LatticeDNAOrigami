// utility.h

#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <array>
#include <random>
#include <functional>

using std::vector;
using std::array;
using std::bind;

namespace Utility {
    
    // Find index for given element
//    template<typename Element_T>
//    int index(vector<Element_T> container, Element_T element);
    int index(vector<int> container, int element);

    // Exception types
    struct NoElement {};
    struct IndexOutOfRange {};
    struct OrigamiMisuse {};
    struct NotImplemented {};
    struct MoveRejection {};

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

        private:
            array<int, 3> m_container;
    };

    // Hash needs this, ambiguous to also have method
    bool operator==(const VectorThree& v1, const VectorThree& v2);

    // Base unit vectors
    const VectorThree xhat {1, 0, 0};
    const VectorThree yhat {0, 1, 0};
    const VectorThree zhat {0, 0, 1};

    // All possible unit vectors
    vector<VectorThree> vectors {
            {1, 0, 0},
            {-1, 0, 0},
            {0, 1, 0},
            {0, -1, 0},
            {0, 0, 1},
            {0, 0, -1}};

    // Random numbers
    std::mt19937_64 random_engine {};
    std::uniform_real_distribution<double> uniform_real_dist {};
    auto gen_uniform_real = bind(uniform_real_dist, random_engine);
    std::uniform_int_distribution<int> uniform_int_dist_vectors(0, 5);
    auto gen_uniform_vector_i = bind(uniform_int_dist_vectors, random_engine);

}


#endif // UTILITY_H
