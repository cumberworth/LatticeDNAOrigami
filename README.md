# Monte Carlo simulation program for a DNA origami lattice model

A program for running Monte Carlo (MC) simulations of a lattice model of DNA origami, especially for studying its self-assembly process.

The original version of the model and simulation techniques implemented here are described originally in [Ref. 1](), while an updated version can be found in [Ref. 2]().
Simulations are run in the grand canonical ensemble, with annealing, replica exchange, and umbrella sampling variants available and configurable through a plain text input parameter file.
Relatively general facilities are available for defining order parameters and associated bias functions in JSON formatted files.
Output is in simple plain text file and JSON formats, and real-time visualization is possible with VMD via the bundled TCL scripts.

This is the first program I wrote in C++, and it has been developed during and after my PhD.
I have tried to follow good practices, but my understanding of what makes good practice has evolved significantly over time, and some design decisions made, particularly regarding the inheritance system for the move type classes, can make modifying and extending the code more challenging than would be optimal.
If anyone is seriously interested in working with this code base, I would be happy to put time in to better modernize and document relevant parts of the code.

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
To run a simulation with `configuration_file.inp`, run
```
latticeDNAOrigami -i [configuration_file.inp]
```
For running a simulation with MPI on `procs` processors, run
```
mpirun -np [procs] latticeDNAOrigami -i [configuration file]
```
When using MPI, it is especially important to check the validity of all configuration file options, as informative exceptions are generally not output.
Consider trying to run a single core version of the desired simulation if it throws MPI errors at the start of the simulation.

## Analysis and visualization

fix this, maybe say something about the types of files that are output

A python package for analyzing the results of simulations (origamipy) is also provided.
Example scripts using the package are located in

`scripts/analysis/`

Configurations can be visualized in two ways.
For vector based graphics, latex scripts using pgf/tikz library are provided in

`scripts/tikz/`

Alternatively, VMD may be used with the scripts provided in

`scripts/vmd/`

## Configuration file specification

fill this out
To see all configuration options, run
```
latticeDNAOrigami -h
```

## System input file specification

An example file is given in `examples/snodin_unbound.json`.
The structure is as follows
```
origami:
    identities: list of lists
    sequences: list of lists
    cyclic: bool
    configurations: list of objects
        step: int
        chains:
            identity: int
            positions: list of lists
            orientations: list of lists
```

### `identities`

This is a list of the identities of each domain in each chain.
The first index in the list is the scaffold, and subsequent lists are the staples.
The scaffold list should contain negative integers, and the staples positive integers.
The index in the list for each corresponds to the domains along that chain.
Domains that are complementary should have identities that sum to 0.
This list should be consistent with the list of sequences, but if using a uniform potential where the energies are set in the configuration file, only the identities are required.

### `sequences`

A list of lists of the sequences of each domain.
Sequences should be given 5' to 3'.

### `cyclic`

Specify if the scaffold is cyclic.

### `configurations`

A list of objects containing a series of configurations for the system.
Usually only a single configuration is included as, which can be used as a starting configuration for the simulation.

#### `step`

The step number, or simply the configuration index.

#### `chains`

A list of objects, where each object is a chain in the system.

##### `identity`

The identity of the chain.
This is the index of the chain in `identities`.

##### `positions`

A list of the coordinates for each domain in the chain.
For each domain, the list is the x, y, and z components of the position vector.

##### `orientations`

A list of the orientations for each domain in the chain.
For each domain, the list is the x, y, and z components of the orientation vector.

## Move set file specification

An example file is given in `examples/moveset_standard.json`.
The structure is as follows
```
origami:
    movetypes: list of ojects
        label: string
        freq: string
        type: string
        [move type specific options]
```

### `label`

A label for the move type; will be used in the output files but is not otherwise used by the program.

### `freq`

A fraction that specifies the frequency of this move type.
It must be a string of a fraction, with no spaces (e.g. `1/2`).

### `type`

The class of move type.
Below are listed the types, along with move type specific options.
See [Ref ?]() for definitions and details of the move types.

#### `OrientationRotation`

Proposes a new orientation vector with uniform probability.

#### `MetStapleExchange`

Proposes either addition or removal of a staple, and uses the symmetric growth scheme.

##### `adaptive_exchange`

Experimental feature; set to false.

#### `MetStapleRegrowth`

Proposes a new configuration for a random staple in the system, and uses the symmetric regrowth scheme to regrow from one of the bound domains.

#### `CBStapleRegrowth`

Proposes a new configuration for a random staple in the system, and uses the configurational bias (CB) regrowth scheme to regrow from one of the bound domains.

#### `CTCBScaffoldRegrowth`

Proposes a new configuration for the scaffold with the conserved topology configurational bias (CTCB) regrowth scheme with contiguous single-segment selection.

##### `max_regrowth`

The maximum number of scaffold domains to regrow.

#### `CTCBJumpScaffoldRegrowth`

Proposes a new configuration for the scaffold with the CTCB regrowth scheme with contiguous multiple-segment selection.
It has the following options in addition to those of `CTCBScaffoldRegrowth`.

##### `max_seg_regrowth`

The maximum number of scaffold domains to regrow in a single segment.

#### `CTRGScaffoldRegrowth`

Proposes a new configuration for the scaffold with the conserved topology recoil growth (CTCB) regrowth scheme with contiguous single-segment selection.
It has the following options in addition to those of `CTCBScaffoldRegrowth`.

##### `max_num_recoils`

The maximum number of recoils allowed.

##### `max_c_attempts`

The maximum number of configurations to try for a given domain.

#### `CTRGJumpScaffoldRegrowth`

Proposes a new configuration for the scaffold with the CTRG scheme with contiguous multiple-segment selection.
It has the sum of the sets of options of `CTCBJumpScaffoldRegrowth` and `CTRGScaffoldRegrowth`.

## Order parameter file specification

All order parameters are functions of the origami system that return single integers.
An example file is given in `examples/ops_standard.json`.
The structure is as follows
```
origami:
    order_params: list of objects
        label: string
        tag: string
        level: int
        ops: list of strings
        type: string
        [order parameter type specific options]
```

### `label`

Label used only for organization of the order parameters file itself.

### `tag`

This is the label that is actually used in the output, and for reference in this or other input files.
It should not contain spaces.

### `level`

Specifies what level in the definition hierarchy the order parameter is.
The number indicates the dependency depth.
For example, level 0 corresponds to a base order parameter that depends on no other order parameter, while level 1 corresponds to an order parameter that depends on a level 0 order parameter.

##### `ops`

A list of order parameter tags that this order parameter directly depends on.
If the order parameter does not depend on any other order parameters, leave empty or do not include.

### `type`

The class of order parameter.
Below are listed the types, along with order parameter type specific options.

#### `NumStaples`

The number staples bound to the system.
This is a base order parameter.

#### `NumStaplesType`

The number of staples bound to the system of type `staple`.
This is a base order parameter.

##### `staple`

The chain index as defined by position in `identities` of the system input file.

#### `StapleTypeFullyBound`

The number of fully bound staples of type `staple`.
It has the same options as `NumStaplesType`.
This is a base order parameter.

#### `NumBoundDomainPairs`

The number of bound domain pairs.
This is a base order parameter.

#### `NumMisboundDomainPairs`

The number of misbound domain pairs.
This is a base order parameter.

#### `NumStackedPairs`

The number of stacked pairs.
This is a base order parameter.

#### `Dist`

The distance, in lattice sites, between two domains.

##### chain1

Index of the chain that the first domain is in.

##### chain2

Index of the chain that the second domain is in.

##### domain1

Index of the first domain in its chain.

##### domain2

Index of the second domain in its chain.

#### `AdjacentSite`

Check for whether two domains are on adjacent sites.
It takes a value of if 1 true, and if 0 false.
It has the same options as `NumStaplesType`.

#### `Sum`

The sum of the ops listed in `ops`.

## Bias function file specification

All bias functions are functions of the origami system that return an energy ().
An example file is given in `examples/biases_mwus-numfulldomains.json`.
The structure is as follows
```
origami:
    bias_functions: list of objects
        label: string
        tag: string
        level: int
        ops: list of strings
        bias_funcs: list of strings
        type: string
        [bias function specific options]  
```

### `label`

Label used only for organization of the bias functions file itself.

### `tag`

This is the label that is actually used in the output, and for reference in this or other input files.
It should not contain spaces.

### `level`

Specifies what level in the definition hierarchy the bias function is.
The number indicates the dependency depth.
For example, level 0 corresponds to a base bias function that depends on no other bias functions, while level 1 corresponds to a bias function that depends on a level 0 bias function.
Note that this is independent from order parameter levels.

### `ops`

A list of order parameter tags that the bias function directly depends on.
If the bias function does not depend on any order parameters, leave empty or do not include.

### `bias_funcs`

A list of bias function tags that this bias function directly depends on.
If the bias function does not depend on any other bias functions, leave empty or do not include.

### `type`

The class of order parameter.
Below are listed the types, along with order parameter type specific options.

#### `LinearStep`

0 at and below `min_op`, `max_bias` at and above `max_op`, and linear between.
Requires one op tag in `ops`.

#### `LinearStepWell`

Equal to `well_bias` between and at `min_op` and `max_op` and linearly increasing outside.
If `op` < `min_op`, then the bias is equal to `slope` * (`min_op` - `op` - 1) + `min_bias`.
If `op` > `max_op`, then the bias is equal to `slope` * (`op` - `max_op` - 1) + `min_bias`.

#### `SquareWell`

If `op` < `min_op` or `op` > `max_op`, then equal to `outside_bias`.
Else, equal to `well_bias`.

#### `Grid`

Biases are a grid with dimensions equal to the number of order parameters listed.
The bias for a given combination of values for the order parameters it is defined in is equal to the value stored in the grid.
This is used by umbrella sampling.

## MWUS window file specification

## Analysis and visualization

## Links
