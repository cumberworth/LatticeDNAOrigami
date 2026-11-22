# Integration tests

`run_tests.py` (Python standard library only) runs `latticeDNAOrigami` on small systems and checks the results. CTest runs it automatically (see `../CMakeLists.txt`).

| Test | What it checks | Catches |
|---|---|---|
| `enumeration` | Exact enumeration of the four-domain system reproduces the stored weights | Changes to the potential or enumeration code |
| `mc-enumeration` | MC order-parameter distributions (bound pairs, stacked pairs, staples, misbound pairs) agree with exact enumeration | Detailed-balance errors, wrong acceptance factors, missing or double-counted energy terms |
| `consistency` | Runs with a full energy and constraint recheck after every step (`constraint_check_freq=1`); the program throws on any drift | Bookkeeping errors in incremental energy updates |
| `regression` | A fixed-seed run reproduces the stored trajectory exactly (order parameters, energies, configurations) | Any behaviour change in refactors that should be behaviour-neutral |

## Running

With CMake presets (see `CMakePresets.json` at the repository root):
```
cmake --preset dev && cmake --build --preset dev
ctest --preset quick       # unit + quick integration tests, ~5 min on 3 cores
ctest --preset slow        # thorough MC checks (use many cores, e.g. the cloud batch)
cmake --preset sanitize && cmake --build --preset sanitize
ctest --preset sanitize    # consistency, regression and enumeration under ASan/UBSan
```
A single test can also be run directly, e.g. `run_tests.py --exe build/dev/latticeDNAOrigami mc-enumeration --temp 345 --seeds 3`.

## Systems and temperatures

- **Four-domain system** (`systems/four_unbound.json`, from `examples/`): two two-domain staples, at most 2 staples.
  - At 340 K it splits between half and fully assembled states (stacked pairs 0–2), which tests the second staple and stacking.
  - At 345 K it splits between empty and half-assembled states, which tests insertion and first binding.
- **System S** (`systems/snodin_unbound.json`, 24 domains): used for the consistency and regression runs, where a larger state space exercises more code paths.

## Precision

The quick MC tests (3 seeds × 7×10⁵ steps) have standard errors of about 2–4% on the probabilities, and fail beyond 5 standard errors. That catches gross errors only. The thorough tests (8 seeds × 10⁷ steps) reach about 1%. Run them before trusting a new or changed move type. To test a new move type, add it to a moveset in `systems/` and pass `--moveset`.

## References

Files in `reference/` were generated on 2026-11-22 from the unmodified code at commit `2ae5905` (the paper's v1.0.0 plus typo fixes). The toolchain was GCC 16.2.1, Boost 1.92 and OpenMPI 5.0.11 on Arch Linux.

- **Enumeration weights** should be portable. They are compared to a relative tolerance of 10⁻⁵.
- **Regression trajectories** depend on the toolchain: floating-point details and libstdc++'s distribution algorithms. On a different compiler or platform, regenerate them from a build known to be correct, using `run_tests.py ... regression --update` with the same arguments as in `CMakeLists.txt`.
- The regression runs were checked to give identical trajectories with the jemalloc allocator, so they don't depend on memory addresses. Changes that only alter allocation patterns (for example replacing containers) can therefore be verified bit for bit.
