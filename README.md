Simulation package for lattice models of the self-assembly of DNA origami. The models are intended to include only the most important details for allowing an explicit representation of the geometry of the system. Simulations are run in the grand cannonical ensemble, with constant temperature, annealing, and temperature parallel tempering variants available. Output is in simple text files, and visulation is accomplished with tikz/pgfplots scripts.

Implementation is with C++11.

Dependencies:

    Libraries:
        boost
        JsonCpp
        openmpi (other backends may be used by the boost mpi wrapper)

    Standalone programs:
        latex

    Latex packages:
        pgf/tikz
        tikz-3d-plot

Running simulations:

An example scripts are given as tests/*.inp. To run a serial simulation:

latticeDNAOrigami -i [simulation script] > [log file]

To run a parallel simulation:

mpirun -np [procs] latticeDNAOrigami -i [simulation script] > [log file]

Visualizing simulations:

Analyzing simulations:
