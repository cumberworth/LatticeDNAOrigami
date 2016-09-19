// utility.h

#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <array>
#include <functional>
#include <unordered_map>
#include <memory>

using std::vector;
using std::array;
using std::bind;
using std::unordered_map;
using std::pair;
using std::unique_ptr;

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

    // Constrants
    const double NA {6.022140857};

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
            int sum();
            int abssum();

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
    const vector<VectorThree> vectors {
            {1, 0, 0},
            {-1, 0, 0},
            {0, 1, 0},
            {0, -1, 0},
            {0, 0, 1},
            {0, 0, -1}};
}

#endif // UTILITY_H
