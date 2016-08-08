// utility.h

#ifndef UTILITY_H
#define UTILITY_H

#include <vector>
#include <array>

using std::vector;
using std::array;

namespace Utility {
    
    // Find index for given element
//    template<typename Element_T>
//    int index(vector<Element_T> container, Element_T element);
    int index(vector<int> container, int element);

    // Exception types
    struct NoElement {};
    struct IndexOutOfRange {};
    struct OrigamiMisuse {};

    // Container for chain domain indices
    struct CDPair {
        int c;
        int d;
        bool operator==(const CDPair& cd_j) {return (this->c == cd_j.c and
               this->d == cd_j.d);};
    };

    inline bool operator==(const CDPair& cd_i, const CDPair& cd_j) {return (cd_i.c == cd_j.c
            and cd_i.d == cd_j.d);};

    // Vector for Z3 space
    class VectorThree: public array<int, 3> {
        public:
            VectorThree(int x, int y, int z): array<int, 3> {{x, y, z}} {};
            VectorThree(): array<int, 3> {{}} {};
            VectorThree operator+(const VectorThree& v_2);
            VectorThree operator-(const VectorThree& v_2);
            VectorThree operator-();
            VectorThree rotate_half(VectorThree axis);
    };

    // Unit vectors
    const VectorThree xhat {1, 0, 0};
    const VectorThree yhat {0, 1, 0};
    const VectorThree zhat {0, 0, 1};
}


#endif // UTILITY_H
