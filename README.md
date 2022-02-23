# Monte Carlo simulation program for a DNA origami lattice model

A program for running Monte Carlo (MC) simulations of a lattice model of DNA origami, especially for studying its self-assembly process.

The original version of the model and simulation techniques implemented here are described originally in [Ref. 1](), while an updated version can be found in [Ref. 2]().
Simulations are run in the grand canonical ensemble, with annealing, replica exchange, and umbrella sampling variants available and configurable through a plain text input parameter file.
Relatively general facilities are available for defining order parameters and associated bias functions in JSON formatted files.
Output is in simple plain text file and JSON formats, and real-time visualization is possible with VMD via the bundled TCL scripts.

## Installation

The program requires the Boost library and an MPI implementation compatible with Boost.MPI.
Building and installation is done with CMake.
To configure, build and install in `installdir`, run
```
CXX=clang++ cmake -S . -B build -DCMAKE_INSTALL_PREFIX=[installdir] -DBUILD_TESTING=OFF
cmake --build build
cmake --install build
```

## Running a simulation

Examples of configuration files and scripts for creating instances of input files can be found in `scripts/simutils/`.
To see all configuration options, run
```
latticeDNAOrigami -h
```
To run a simulation with `configuration_file.inp`, run
```
latticeDNAOrigami -i [configuration_file.inp]
```
For running a simulation with MPI on `procs` processors, run
```
mpirun -np [procs] latticeDNAOrigami -i [configuration file]
```

## Analysis and visualization

A python package for analyzing the results of simulations (origamipy) is also provided. Example scripts using the package are located in

`scripts/analysis/`

Configurations can be visualized in two ways.
For vector based graphics, latex scripts using pgf/tikz library are provided in

`scripts/tikz/`

Alternatively, VMD may be used with the scripts provided in

`scripts/vmd/`

## Links
