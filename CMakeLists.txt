cmake_minimum_required(VERSION 3.15)

project(
  LatticeDNAOrigami
  VERSION 0.1.0
  LANGUAGES CXX)

# Build type
set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${DEFAULT_BUILD_TYPE}"
      CACHE STRING "Choose the type of build." FORCE)

  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(cmake/CheckGit.cmake)
checkgitsetup()

add_library(
  LatticeDNAOrigami_lib
  src/annealing_simulation.cpp
  src/bias_functions.cpp
  src/cb_movetypes.cpp
  src/constant_temp_simulation.cpp
  src/domain.cpp
  src/enumerate.cpp
  src/files.cpp
  src/hash.cpp
  src/ideal_random_walk.cpp
  src/jsoncpp.cpp
  src/met_movetypes.cpp
  src/movetypes.cpp
  src/nearest_neighbour.cpp
  src/order_params.cpp
  src/orientation_movetype.cpp
  src/origami_potential.cpp
  src/origami_system.cpp
  src/parser.cpp
  src/ptmc_simulation.cpp
  src/random_gens.cpp
  src/rg_movetypes.cpp
  src/simulation.cpp
  src/top_constraint_points.cpp
  src/transform_movetypes.cpp
  src/us_simulation.cpp
  src/utility.cpp)

target_compile_features(LatticeDNAOrigami_lib PRIVATE cxx_std_17)
target_include_directories(LatticeDNAOrigami_lib PUBLIC include)

target_link_libraries(LatticeDNAOrigami_lib PUBLIC git_version)

# MPI
find_package(MPI REQUIRED)
message(
  STATUS
    "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS"
)
target_link_libraries(LatticeDNAOrigami_lib PUBLIC MPI::MPI_CXX)

# Boost
find_package(Boost 1.69 REQUIRED COMPONENTS program_options mpi serialization
                                            system filesystem)
message(STATUS "Boost version: ${Boost_VERSION}")
target_link_libraries(
  LatticeDNAOrigami_lib
  PUBLIC Boost::filesystem Boost::program_options Boost::mpi
         Boost::serialization Boost::system Boost::filesystem)

# Interprocedular optimization
include(CheckIPOSupported)
check_ipo_supported(RESULT RESULT)
if(RESULT)
  message(STATUS "Using interprocedural optimization.")
  set_target_properties(LatticeDNAOrigami_lib
                        PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# This used to be needed

# Brute force method is required to add interprocedural optimization because
# CMake does not support it for clang There is a bug in the llvm code that
# prevents this from being used in combination with clang-tidy
# target_compile_options(LatticeDNAOrigami_lib PUBLIC -flto=thin)
# set_target_properties(LatticeDNAOrigami_lib PROPERTIES LINK_FLAGS
# "-fuse-ld=lld -flto=thin")

# CCache (increased compiliation speed)
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# Compiler warnings

# target_compile_options(LatticeDNAOrigami_lib PUBLIC -Wall -Wextra -Wshadow
# -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
# -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference
# -Wdouble-promotion -Werror)

# target_compile_options(LatticeDNAOrigami_lib PUBLIC -Wall -Wextra -Wshadow
# -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Woverloaded-virtual
# -Wpedantic -Wconversion -Wsign-conversion -Wnull-dereference
# -Wdouble-promotion)

# Create program
add_executable(latticeDNAOrigami apps/main.cpp)
target_link_libraries(latticeDNAOrigami PUBLIC LatticeDNAOrigami_lib)

# Testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(tests)
endif()

# Install executable
install(TARGETS latticeDNAOrigami RUNTIME DESTINATION bin)
