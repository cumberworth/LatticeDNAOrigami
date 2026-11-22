#!/usr/bin/env python3
"""Integration tests for latticeDNAOrigami (standard library only).

Tests:
  enumeration      exact enumeration reproduces the stored reference weights
  mc-enumeration   Monte Carlo order-parameter distributions agree with exact
                   enumeration (catches detailed-balance and energy errors)
  consistency      the full energy/constraint recheck passes on every check
                   during a run (the program exits non-zero otherwise)
  regression       a fixed-seed run reproduces the stored reference trajectory
                   exactly (for refactors that should not change behaviour)

Usage:
  run_tests.py --exe PATH/latticeDNAOrigami <test> [options]
  run_tests.py --exe ... regression --update      # regenerate a reference

References in reference/ were generated with the build described in
reference/README.md. The regression references depend on the toolchain
(floating-point details and libstdc++'s distribution algorithms), so regenerate
them with --update when switching compiler or platform, from a build known to be
correct.
"""

import argparse
import math
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SYSTEMS = HERE / "systems"
REFERENCE = HERE / "reference"
OP_TAGS = ["numfulldomains", "nummisdomains", "numstackedpairs", "numstaples"]

BASE_PARAMS = {
    "domain_type": "HalfTurn",
    "binding_pot": "FourBody",
    "misbinding_pot": "Opposing",
    "stacking_pot": "Constant",
    "hybridization_pot": "NearestNeighbour",
    "apply_mean_field_cor": "false",
    "staple_M": "1e-7",
    "cation_M": "0.5",
    "staple_u_mult": "1",
    "stacking_ene": "-1000",
    "domain_update_biases_present": "false",
    "order_parameter_file": "ops_standard.json",
    "ops_to_output": " ".join(OP_TAGS),
}

FOUR = {
    "origami_input_filename": "four_unbound.json",
    "min_total_staples": "0",
    "max_total_staples": "2",
    "max_type_staples": "2",
    "max_staple_size": "2",
}

SNODIN = {
    "origami_input_filename": "snodin_unbound.json",
    "max_total_staples": "24",
    "max_type_staples": "12",
    "max_staple_size": "2",
}


class TestFailure(Exception):
    pass


def write_inp(path, params):
    path.write_text("".join(f"{k}={v}\n" for k, v in params.items()))


def start_program(exe, workdir, params, name):
    """Start latticeDNAOrigami in workdir with the given parameters."""
    for f in SYSTEMS.iterdir():
        shutil.copy(f, workdir)
    (workdir / "outs").mkdir(exist_ok=True)
    params = {**params, "output_filebase": f"outs/{name}"}
    inp = workdir / f"{name}.inp"
    write_inp(inp, params)
    proc = subprocess.Popen([str(exe), "-i", inp.name], cwd=workdir, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)

    return proc, workdir / "outs" / name, name


def finish_program(proc, name):
    output, _ = proc.communicate()
    if proc.returncode != 0:
        tail = "\n".join(output.strip().split("\n")[-15:])
        raise TestFailure(f"{name}: program exited with {proc.returncode}\n{tail}")


def run_program(exe, workdir, params, name):
    proc, out, name = start_program(exe, workdir, params, name)
    finish_program(proc, name)

    return out


def mc_params(temp, steps, seed, check_freq, output_freq, moveset):
    return {
        **BASE_PARAMS,
        "temp": str(temp),
        "simulation_type": "constant_temp",
        "movetype_file": moveset,
        "random_seed": str(seed),
        "ct_steps": str(steps),
        "max_duration": "100000",
        "centering_freq": str(output_freq),
        "constraint_check_freq": str(check_freq),
        "logging_freq": str(max(output_freq, steps // 10)),
        "configs_output_freq": str(output_freq),
        "vtf_output_freq": "0",
        "counts_output_freq": str(output_freq),
        "order_params_output_freq": str(output_freq),
        "times_output_freq": "0",
        "energies_output_freq": str(output_freq),
    }


def read_weights(path):
    weights = {}
    for line in Path(path).read_text().split("\n")[1:]:
        m = re.match(r"\((.*)\) (\S+)", line)
        if m:
            weights[tuple(int(v) for v in m.group(1).split())] = float(m.group(2))

    return weights


def read_ops(path):
    lines = Path(path).read_text().strip().split("\n")
    header = lines[0].replace(",", " ").split()
    rows = [tuple(int(v) for v in line.split()[1:]) for line in lines[1:]]

    return header, rows


def marginal(weights, index):
    p = defaultdict(float)
    for state, w in weights.items():
        p[state[index]] += w

    return p


def batch_means_error(indicator, nbatches=20):
    """Standard error of the mean of a correlated 0/1 series via batch means."""
    n = len(indicator) // nbatches * nbatches
    if n == 0:
        return math.inf
    size = n // nbatches
    means = [sum(indicator[i * size:(i + 1) * size]) / size for i in range(nbatches)]
    mean = sum(means) / nbatches
    var = sum((m - mean) ** 2 for m in means) / (nbatches - 1)

    return math.sqrt(var / nbatches)


# ---------------------------------------------------------------------------
# Tests


def test_enumeration(args, workdir):
    ref = REFERENCE / f"four-{args.temp}K.weights"
    params = {**BASE_PARAMS, **FOUR, "temp": str(args.temp), "simulation_type": "enumerate",
              "enumerate_staples_only": "false"}
    out = run_program(args.exe, workdir, params, f"enum-{args.temp}")
    weights = Path(f"{out}.weights")
    if args.update:
        shutil.copy(weights, ref)
        return f"updated {ref.name}"

    got, expected = read_weights(weights), read_weights(ref)
    if set(got) != set(expected):
        raise TestFailure(f"state sets differ: {sorted(set(got) ^ set(expected))}")
    worst = max(abs(got[s] - expected[s]) / max(expected[s], 1e-300) for s in expected if expected[s] > 1e-12)
    if worst > 1e-5:
        raise TestFailure(f"weights differ from reference (max relative difference {worst:.2e})")

    return f"{len(got)} states match the reference (max relative difference {worst:.1e})"


def test_mc_enumeration(args, workdir):
    """Compare MC marginals of each order parameter with exact enumeration.

    With --seeds 1 the error is estimated by batch means within the run. With
    --seeds N > 1, N independent runs go in parallel and the error comes from
    the spread between them, which is more robust for slowly mixing states.
    Each order-parameter value with exact probability above `min_p` must lie
    within `sigmas` standard errors (with a small floor on the error so
    near-certain states don't produce spurious failures). Seeds are fixed, so
    the outcome is deterministic for a given build.
    """
    ref = read_weights(REFERENCE / f"four-{args.temp}K.weights")
    runs = []
    procs = []
    for i in range(args.seeds):
        run_dir = workdir / f"seed-{i}"
        run_dir.mkdir()
        params = {**BASE_PARAMS, **FOUR, **mc_params(args.temp, args.steps, args.seed + i, 10000,
                                                       args.output_freq, args.moveset)}
        procs.append(start_program(args.exe, run_dir, params, f"mc-{args.temp}-{i}"))
    for proc, out, name in procs:
        finish_program(proc, name)
        header, rows = read_ops(f"{out}.ops")
        runs.append((header, rows[len(rows) // 10:]))  # discard burn-in

    report, failures = [], []
    for tag in ("numfulldomains", "numstackedpairs", "numstaples", "nummisdomains"):
        i_ref = OP_TAGS.index(tag)
        exact = marginal(ref, i_ref)
        for value, p in sorted(exact.items()):
            if p < args.min_p:
                continue
            per_run = []
            for header, rows in runs:
                i_mc = header.index(tag)
                per_run.append([1 if r[i_mc] == value else 0 for r in rows])
            pooled = [x for ind in per_run for x in ind]
            p_mc = sum(pooled) / len(pooled)
            if len(per_run) > 1:
                means = [sum(ind) / len(ind) for ind in per_run]
                var = sum((m - p_mc) ** 2 for m in means) / (len(means) - 1)
                err = math.sqrt(var / len(means))
            else:
                err = batch_means_error(pooled)
            err = max(err, 0.003)
            z = (p_mc - p) / err
            line = f"{tag}={value}: exact {p:.4f}, MC {p_mc:.4f} ± {err:.4f} (z = {z:+.1f})"
            report.append(line)
            if abs(z) > args.sigmas:
                failures.append(line)

    if failures:
        raise TestFailure("MC disagrees with enumeration:\n  " + "\n  ".join(failures) + "\nall:\n  "
                          + "\n  ".join(report))

    nsamples = sum(len(rows) for _, rows in runs)
    return f"{nsamples} samples from {args.seeds} run(s) agree with enumeration\n  " + "\n  ".join(report)


def test_consistency(args, workdir):
    system = FOUR if args.system == "four" else SNODIN
    params = {**BASE_PARAMS, **system, **mc_params(args.temp, args.steps, args.seed, 1, args.output_freq,
                                                     args.moveset)}
    run_program(args.exe, workdir, params, f"consistency-{args.system}")

    return f"{args.steps} steps with a full energy/constraint recheck after every step"


def test_regression(args, workdir):
    system = FOUR if args.system == "four" else SNODIN
    name = f"regression-{args.system}-{args.temp}K"
    params = {**BASE_PARAMS, **system, **mc_params(args.temp, args.steps, args.seed, 1000, args.output_freq,
                                                     args.moveset)}
    out = run_program(args.exe, workdir, params, name)
    produced = {ext: Path(f"{out}.{ext}") for ext in ("ops", "ene", "trj")}
    if args.update:
        for ext, path in produced.items():
            shutil.copy(path, REFERENCE / f"{name}.{ext}")
        return f"updated {name}.*"

    for ext, path in produced.items():
        expected = (REFERENCE / f"{name}.{ext}").read_text().split("\n")
        got = path.read_text().split("\n")
        for i, (a, b) in enumerate(zip(expected, got)):
            if a != b:
                raise TestFailure(f"{name}.{ext} first differs at line {i + 1}:\n  expected: {a[:120]}\n  got:      {b[:120]}")
        if len(expected) != len(got):
            raise TestFailure(f"{name}.{ext}: {len(got)} lines, expected {len(expected)}")

    return f"{name}: ops, energies and configurations identical to the reference"


TESTS = {
    "enumeration": test_enumeration,
    "mc-enumeration": test_mc_enumeration,
    "consistency": test_consistency,
    "regression": test_regression,
}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--exe", required=True, type=Path, help="Path to latticeDNAOrigami")
    parser.add_argument("test", choices=TESTS)
    parser.add_argument("--temp", type=int, default=340)
    parser.add_argument("--system", choices=["four", "snodin"], default="four")
    parser.add_argument("--moveset", default="moveset_four.json")
    parser.add_argument("--steps", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=20261122)
    parser.add_argument("--seeds", type=int, default=1, help="Independent runs in parallel (mc-enumeration)")
    parser.add_argument("--output_freq", type=int, default=100)
    parser.add_argument("--min_p", type=float, default=0.02, help="Smallest exact probability tested")
    parser.add_argument("--sigmas", type=float, default=5.0, help="Allowed deviation in standard errors")
    parser.add_argument("--update", action="store_true", help="Regenerate the stored reference")
    parser.add_argument("--keep", action="store_true", help="Keep the working directory")
    args = parser.parse_args()
    args.exe = args.exe.resolve()

    workdir = Path(tempfile.mkdtemp(prefix=f"ldo-{args.test}-"))
    try:
        message = TESTS[args.test](args, workdir)
        print(f"PASS {args.test}: {message}")
        status = 0
    except TestFailure as e:
        print(f"FAIL {args.test}: {e}")
        status = 1
    finally:
        if args.keep:
            print(f"(working directory kept: {workdir})")
        else:
            shutil.rmtree(workdir, ignore_errors=True)

    sys.exit(status)


if __name__ == "__main__":
    main()
