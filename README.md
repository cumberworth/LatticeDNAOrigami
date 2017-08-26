Simulation package for lattice models of the self-assembly of DNA origami. The models are intended to include only the most important details for allowing an explicit representation of the geometry of the system. Simulations are run in the grand canonical ensemble, with annealing, Hamiltonian parallel tempering, and umbrella sampling variants available and configurable through a plain text input parameter file. Relatively general facilities are available for defining order parameters and associated bias functions in JSON formatted files. Output is in simple plain text file and JSON formats, and real-time visualization is possible with VMD via the bundled TCL scripts.

The core simulation package is implemented with C++14.

Core dependencies:

* boost
* JsonCpp
* openmpi (other backends may be used by the boost mpi wrapper)

Running simulations:

Example scripts can be found in `scripts/simulation_examples/*.inp`. To see all configuration options, run

`latticeDNAOrigami -h`

 To run a simulation, enter

latticeDNAOrigami -i [configuration file] > [log file]

To run a parallel simulation, enter

mpirun -np [procs] latticeDNAOrigami -i [configuration file] > [log file]
