find_package(Catch2 REQUIRED)
include(CTest)
include(Catch)

add_executable(tests src/test_main.cpp src/test_nearest_neighbour.cpp
                     src/test_ideal_random_walk.cpp)
target_include_directories(tests PUBLIC ../include)
target_link_libraries(tests Catch2::Catch2)
target_link_libraries(tests LatticeDNAOrigami_lib)

catch_discover_tests(tests)
