#include <catch.hpp>

#define private public
#define protected public

#include "ideal_random_walk.h"

using std::cout;

using namespace Utility;
using namespace IdealRandomWalk;

SCENARIO("Computed ideal random walks match hand calculated values") {
    IdealRandomWalks ideal_walk {};
    VectorThree origin {0, 0, 0};
    VectorThree start_pos;
    VectorThree end_pos;
    int N;
    long double number_of_walks;

    GIVEN("A typical difference in positions and number of steps") {
        start_pos = origin;
        end_pos = {0, 2, 3};
        N = 5;
        number_of_walks = 10;
        REQUIRE(number_of_walks == ideal_walk.num_walks(start_pos, end_pos, N));

        start_pos = origin;
        end_pos = {0, 2, 3};
        N = 7;
        number_of_walks = 665;
        REQUIRE(number_of_walks == ideal_walk.num_walks(start_pos, end_pos, N));
    }

    GIVEN("No remaining steps") {
        start_pos = origin;
        end_pos = {0, 2, 3};
        N = 6;
        number_of_walks = 0;
        REQUIRE(number_of_walks == ideal_walk.num_walks(start_pos, end_pos, N));
    }

    GIVEN("A typical difference in positions and large number of steps") {
        start_pos = origin;
        end_pos = {0, 2, 3};
        N = 51;
        number_of_walks = 5.947398897268465e+36;
        REQUIRE(Approx(number_of_walks) == ideal_walk.num_walks(start_pos, end_pos, N));
    }

    GIVEN("Walks with the same absolulte value of the position difference") {
        start_pos = {1, 2, 3};
        end_pos = {3, 4, 5};
        N = 8;
        long double num_walks_1 {ideal_walk.num_walks(start_pos, end_pos, N)};
        long double num_walks_2 {ideal_walk.num_walks(end_pos, start_pos, N)};
        VectorThree DR {2, 2, 2};
        long double num_walks_3 {ideal_walk.num_walks(origin, DR, N)};

        THEN("The number of walks is the same") {
            REQUIRE(num_walks_1 == num_walks_2);
            REQUIRE(num_walks_1 == num_walks_3);
        }
    }

    GIVEN("All permutations of signs on the position difference coordinates") {
        vector<long double> list_num_walks {};
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, 2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {-1, 2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, -2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, 2, -3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {-1, -2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {-1, 2, -3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, -2, -3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {-1, -2, -3}, 8));

        THEN("The number of walks is the same") {
            for (auto n_walks: list_num_walks) {
                REQUIRE(n_walks == list_num_walks[0]);
            }
        }
    }

    GIVEN("All permutations of a set of position difference coordinates") {
        vector<long double> list_num_walks {};
        list_num_walks.clear();
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, 2, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {2, 1, 3}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {3, 2, 1}, 8));
        list_num_walks.push_back(ideal_walk.num_walks(origin, {1, 3, 2}, 8));

        THEN("The number of walks is the same") {
            for (auto n_walks: list_num_walks) {
                REQUIRE(n_walks == list_num_walks[0]);
            }
        }
    }
}
