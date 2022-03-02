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

The configuration file is a simple key-value format.
To see all configuration options, run
```
latticeDNAOrigami -h
```
Example configuration files can be found in `examples/`.
There is an example of a constant temperature simulation, `constant-temp.inp`, a temperature replica exchange simulation, `ptmc.inp`, (also known as temperature parallel tempering, hence the abbreviation) a replica-exchange multi-window umbrella sampling simulation, and an enumeration run, `enum.inp`.
Some additional details of the key parameters are given below

### domain_type

This can be either `HalfTurn` or `ThreeQuarterTurn`.
`HalfTurn` is what is used in [Ref. 1](), [Ref. 2](), and?.
It is used to model systems that have some number of turns plus an additional half turn.
For example, 16 base pairs (bp) corresponds to about 1.5 turns of a helix.
`ThreeQuarterTurn` can be used to model systems with binding domains that correspond to some number of turns plus an additional three-quarter turn of a helix.
For example, 8 bp corresponds to about 0.75 of a helix.
Note that single binding domains in the design can be simulation with multiple binding domains as defined in the model.
For example, 16 bp binding domains can be modeled as a single `Halfturn` binding domain, or two `ThreeQuarterTurn` binding domains.
The `ThreeQuarterTurn` binding domain type also allows to model designs that are 3D (note that the model is always in 3D, but `HalfTurn` binding domains can only produce planar designs).

### binding_pot

The potential used for describing the binding of two complementary binding domains.
Only one option is provided, `FourBody`, which is the potential used and described in [Ref. 2]().
The potential used in [Ref. 1]() is no longer provided, but the code used for that paper can be found in [Ref. 4](), as well as in the history of this repository.

### misbinding_pot

The potential used for describing the misbinding of two binding domains that are not fully complementary.
The potential used in Ref. 1, Ref. 2, and Ref. 3 is `Opposing`, but to prevent misbinding altogether, use `Disallowed`.

### stacking_pot

The potential used to describe helical stacking of adjacent binding domains.
This should be set to `Constant`; `SequenceSpecific` has not been fully implemented.
This means that `stacking_ene` must be set.

### hybridization_pot

The potential used to calculate the hybridization free energy.
Select `NearestNeighbour` to use the nearest neighbour model, which requires the system definition file include sequences for all binding domains.
For a single enthalpy and entropy for all bound domains, and another set for all misbound domains, select `Uniform`.
In this case, `binding_h`, `misbinding_h`, `binding_s`, and `misbinding_s` must be specified.

### `max_total_staples`, `max_type_staples`, `max_staple_size`

These should be set to smallest values that are reasonable.
This is because some of the output file formats require a constant number of staples, and so if these values are unnecessarily big, so too will the output files be.

### `domain_update_biases_present`

This should be set to true if any bias functions being used are defined for specific binding domains in the system.
For example, if a bias function is defined on the distance between two binding domains.
This ensures that the bias is updated each time a domain is set, rather than after the whole configuration of the system has been set.

### `simulation_type`

The type of simulation to run.

#### `enumerate`

This will enumerate all possible configurations for the system for the given limits on the numbers of staples.
It is only possible for relatively small systems (i.e., scaffolds up to 5 binding domains).
With the option `enumerate_staples_only`, only the configurations of the staples will be enumerated, which also means larger systems are feasible.
This mode is used to validate the move types obey detailed balance by comparing the exact enumeration results to results of simulations on the same system.

#### `constant_temp`

These are standard MC simulations at a constant temperature.
The simulations can be limited by the total number of steps with `ct_steps`, or by the duration in seconds with `max_duration`.

#### `annealing`

#### `[t, ut, hut, st, 2D]_parallal_tempering`

All of `temps`, `chem_pot_mults`, `bias_mults`, and `stacking_mults` must be set, with the number of entries in each (space delimited) equal to the number of replicas `num_reps`.
The temperatures are set with a space delimited list for `temps`.
Multipliers to modify the chemical potential for each replica can be set with `chem_pots_mults`.
Multipliers to modify the total bias for each replica can be set with `bias_mults`.
Multipliers to modify the stacking energy for each replica can be set with `stacking_mults`.
Set the multipliers to an array of ones as a default.
The number of steps between swaps is set with `exchange_interval`.
The simulations can be limited by the number of swaps `swaps` or the duration in seconds `max_duration`.

#### `t_parallel_tempering`

Replica exchange with temperature as the exchange variable.
This means that chemical potential does not change when the temperature does, which means that staple concentration is not constant.

#### `ut_parallel_tempering`

Replica exchange with temperature and chemical potential as a combined exchange variable.
Also exchanging the chemical potential allows the staple concentrations to remain constant.

#### `hut_parallel_tempering`

Replica exchange with temperature, chemical potential, and the bias multiplier as a combined exchange variable.

#### `st_parallel_tempering`

Replica exchange with temperature and the stacking multiplier as a combined exchange variable.

#### `2d_parallel_tempering`

2D replica exchange, where the temperature and the stacking multiplier are the two independent exchange variables.

#### `umbrella_sampling`

#### `mw_umbrella_sampling`

#### `ptmw_umbrella_sampling`

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

This is a list of the identities of each binding domain in each chain.
The first index in the list is the scaffold, and subsequent lists are the staples.
The scaffold list should contain negative integers, and the staples positive integers.
The index in the list for each corresponds to the binding domains along that chain.
Binding domains that are complementary should have identities that sum to 0.
This list should be consistent with the list of sequences, but if using a uniform potential where the energies are set in the configuration file, only the identities are required.

### `sequences`

A list of lists of the sequences of each binding domain.
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

A list of the coordinates for each binding domain in the chain.
For each binding domain, the list is the x, y, and z components of the position vector.

##### `orientations`

A list of the orientations for each binding domain in the chain.
For each binding domain, the list is the x, y, and z components of the orientation vector.

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

The maximum number of scaffold binding domains to regrow.

#### `CTCBJumpScaffoldRegrowth`

Proposes a new configuration for the scaffold with the CTCB regrowth scheme with contiguous multiple-segment selection.
It has the following options in addition to those of `CTCBScaffoldRegrowth`.

##### `max_seg_regrowth`

The maximum number of scaffold binding domains to regrow in a single segment.

#### `CTRGScaffoldRegrowth`

Proposes a new configuration for the scaffold with the conserved topology recoil growth (CTCB) regrowth scheme with contiguous single-segment selection.
It has the following options in addition to those of `CTCBScaffoldRegrowth`.

##### `max_num_recoils`

The maximum number of recoils allowed.

##### `max_c_attempts`

The maximum number of configurations to try for a given binding domain.

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

The distance, in lattice sites, between two binding domains.

##### chain1

Index of the chain that the first binding domain is in.

##### chain2

Index of the chain that the second binding domain is in.

##### domain1

Index of the first binding domain in its chain.

##### domain2

Index of the second binding domain in its chain.

#### `AdjacentSite`

Check for whether two binding domains are on adjacent sites.
It takes a value of if 1 true, and if 0 false.
It has the same options as `NumStaplesType`.

#### `Sum`

The sum of the ops listed in `ops`.

## Bias function file specification

All bias functions are functions of the origami system that return an energy in units of kb.
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

To run MWUS, an additional file is required to specify the windows.
This file is a simple text file.
The first line must contain the tags of the bias functions that will be used to create the restraints for the windows.
These is usually a `LinearStepWell` bias function type.
For 1D MWUS, only one tag is needed.
Each subsequent row defines the minimum and maximum values of the order parameters, which override the `min_op` and `max_op` specified in the definition of the bias function.
If multidimensional MWUS is to be performed, all minimum values are listed first, delimited by spaces, followed by a comma, and then all maximum values, with the order parameters in the same order as the tags in the top row.
An example for 1D MWUS is given in `examples/snodin-numfulldomains.windows`.

## Analysis and visualization

## Links
