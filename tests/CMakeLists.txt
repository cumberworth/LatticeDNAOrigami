find_package(Catch2 3 REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
include(CTest)
include(Catch)

# Unit tests
add_executable(tests src/test_nearest_neighbour.cpp src/test_ideal_random_walk.cpp)
target_include_directories(tests PUBLIC ../include)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain LatticeDNAOrigami_lib)
catch_discover_tests(tests PROPERTIES LABELS "unit")

# Integration tests (see integration/run_tests.py). Labels:
#   quick     run by default; a few minutes in total with ctest -j
#   slow      thorough MC-vs-enumeration checks (~30 min wall each on 8 cores);
#             run with `ctest -L slow`, e.g. before trusting a new move type
set(RUN_TESTS ${CMAKE_CURRENT_SOURCE_DIR}/integration/run_tests.py)

function(add_integration_test name labels)
    add_test(NAME ${name}
             COMMAND ${Python3_EXECUTABLE} ${RUN_TESTS} --exe $<TARGET_FILE:latticeDNAOrigami>
                     ${ARGN})
    set_tests_properties(${name} PROPERTIES LABELS "${labels}")
endfunction()

foreach(temp 340 345)
    add_integration_test(enumeration-${temp}K "quick" enumeration --temp ${temp})
endforeach()

# Quick MC checks: catch gross errors (wrong acceptance factors, missing terms)
foreach(temp 340 345)
    add_integration_test(mc-enumeration-${temp}K "quick"
                         mc-enumeration --temp ${temp} --seeds 3 --steps 700000 --output_freq 1000)
endforeach()
set_tests_properties(mc-enumeration-340K mc-enumeration-345K PROPERTIES PROCESSORS 3)

# Thorough MC checks: ~1% precision on the order-parameter distributions
foreach(temp 340 345)
    add_integration_test(mc-enumeration-thorough-${temp}K "slow"
                         mc-enumeration --temp ${temp} --seeds 8 --steps 10000000 --output_freq 1000)
    set_tests_properties(mc-enumeration-thorough-${temp}K PROPERTIES PROCESSORS 8)
endforeach()

add_integration_test(consistency-four "quick"
                     consistency --system four --temp 345 --steps 200000 --output_freq 1000)
add_integration_test(consistency-snodin "quick"
                     consistency --system snodin --temp 330 --moveset moveset_standard.json
                     --steps 30000 --output_freq 1000)

add_integration_test(regression-four "quick"
                     regression --system four --temp 345 --steps 200000 --output_freq 1000)
add_integration_test(regression-snodin "quick"
                     regression --system snodin --temp 330 --moveset moveset_standard.json
                     --steps 100000 --output_freq 1000)
