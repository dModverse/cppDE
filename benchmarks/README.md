# cppDE benchmark suite

A reproducible comparison of cppDE's own solvers against **SUNDIALS
CVODE(S)**, on published systems-biology models and on the classical
stiff IVP test problems, without sensitivities and with first-order
forward sensitivities.

Nothing here is part of the R package: `benchmarks/` is listed in
`.Rbuildignore` and is never installed, checked, or shipped.

---

## Quick start

```sh
Rscript benchmarks/fetch-models.R                    # once: download the collection
Rscript benchmarks/run-benchmarks.R --tier tiny      # ~5 min
Rscript benchmarks/run-benchmarks.R                  # medium, the default
Rscript benchmarks/run-benchmarks.R --tier full      # everything
```

Results land in `benchmarks/results/<timestamp>/`:

| file | contents |
|---|---|
| `results.csv` | one row per (problem, condition, solver, mode, tolerance) |
| `run-info.txt` | R / cppDE / compiler / CPU, the options used, what was skipped |
| `01…06-*.png` | the figures described below |

`Rscript benchmarks/run-benchmarks.R --help` lists every option.

---

## Three depths, and why the smallest one is trustworthy

| tier | contents | roughly |
|---|---|---|
| `tiny` | one representative of **every** problem trait, 2 tolerances | 5 min |
| `medium` | at least two per trait, wider size and sensitivity range | 30 min |
| `full` | everything within `--max-states` | hours |

`tiny` exists to answer *"did this change help or hurt?"* without waiting
for a full sweep. That only works if it leaves no blind spot: a tier
containing no sparse system cannot see a regression in the KLU path, and
one containing no event cannot see a regression in event handling. So
the problems are not picked by size but by **trait coverage**:

| trait | why it is its own class |
|---|---|
| `stiff-extreme` | rate constants spanning many decades |
| `stiff-moderate` | stiff without extreme scale separation |
| `oscillatory` | limit cycle: error shows up as phase drift |
| `relaxation` | stiffness switches on and off along the trajectory |
| `sparse` | sparse/banded Jacobian: the KLU path |
| `large` | ≥ 100 states |
| `events` | timed discontinuities: stop, jump, restart |
| `many-sens` / `few-sens` | wide versus narrow AD |
| `rational` | Michaelis–Menten / Hill (division by states) |
| `transcendental` | `exp` / `log` / trigonometric terms |
| `log-horizon` | output grid over ≥ 6 decades in time |
| `systems-biology` | published model, realistic parameter values |

Traits are inferred automatically for PEtab models (size, events,
sensitivity count, Jacobian density, division by states, transcendental
calls); the regime traits that cannot be read off the equations
(stiffness, oscillation) are declared on the classical problems.

Tier membership is an **explicit list**, not a computed cover, so that
numbers stay comparable across runs as problems are added. What *is*
computed is the check: every run prints its coverage and warns loudly if
a trait has dropped to zero.

```
Trait coverage [tiny]:
  stiff-extreme       1  scale separation over many decades
  sparse              1  sparse or banded Jacobian -- exercises the KLU path
  events           !! 0  timed discontinuities -- stop, jump, restart
```

`--tol`, `--max-states` and `--max-sens` override the tier's defaults
when passed explicitly.

---

## What is being compared

Both backends are driven through exactly the same pipeline, so the
measurement isolates the integrator and nothing else:

* the **same C++ right-hand side**, generated from the same expression
  strings by the same SymPy codegen;
* the **same analytic Jacobian** (never a finite-difference one);
* the **same output grid, tolerances and initial conditions**;
* sensitivities from **forward AD** (cppDE) versus **CVODES forward
  sensitivity analysis**, in both cases the analytic route, never
  finite differences;
* no R callback on either side.

| solver | what it is |
|---|---|
| `cppDE_ndf` | cppDE's multistep solver, NDF formula (the default) |
| `CVODE_bdf` | SUNDIALS CVODE / CVODES, BDF: the reference implementation |
| `cppDE_bdf` | cppDE with `useNDF = FALSE`, i.e. CVODE's own BDF formula, so a difference isolates the NDF change (`--extra-solvers`) |
| `cppDE_rb4` | Rosenbrock4, a one-step alternative (`--extra-solvers`) |
| `cppDE_dense` / `cppDE_sparse` / `CVODE_dense` / `CVODE_sparse` | the same two backends with the linear solver **pinned** instead of auto-detected, the sparse sweep, see below |

### The sparse sweep is part of the run

Both backends auto-detect sparsity from the symbolic Jacobian
(`decide_sparse()` in `inst/python/codegen_cppODE.py`: sparse from 8
states unless the structure is denser than 0.4). Auto-detection is what
the four solvers above are there to check: every model where the
decision is a real one is additionally run with the linear solver pinned
dense and pinned sparse, on both backends.

This used to be `--sparse-sweep`, a *separate* run that replaced the
head-to-head and filtered the problem list down to the sparse models.
Two runs meant two result folders that could not be read against each
other, the pinned times came from a different process, a different
build, and often a different machine than the auto-detected ones. Now
both live in one `results.csv`, distinguished by the **`pinned`** column
(`auto` / `dense` / `sparse`), and every pinned cell shares its worker
and its neighbourhood with the auto cell it has to be compared to.

| | |
|---|---|
| `--no-sparse-sweep` | leave it out; the cheap run |
| `--max-density <x>` | only sweep models structurally sparser than this (default `0.25`) |
| `--min-sweep-states <n>` | smallest system worth pinning (default `8`, matching the codegen's own threshold) |

The default cutoff sits below the codegen's 0.4, so the sweep measures
*what sparsity buys where it is used*. Raising `--max-density` past 0.4
is what turns it into a test of *whether 0.4 is the right place to
switch*, the two questions want different runs.

The `pinned` column keeps the two apart in every summary: the headline
speed-up reads `auto` rows only, so the sweep can never make a model
count three times in a mean meant to weight it once, and the dense/sparse
table reads pinned rows only, so it never pairs a model against itself.

### Second-order sensitivities

`--modes nosens,sens1,sens2` adds Hessian runtimes. **CVODES has no
second-order sensitivities**, so there is no second implementation to
compare against and none is invented: the `sens2` rows report wall-clock
time and step counts only, with `err` and `err_sens` left `NA`. They are
cppDE measurements, not a comparison.

Second-order forward AD carries M(M+1)/2 additional directions, so cost
grows quadratically in the number of parameters. `--max-sens2` (default
10) caps M for this mode separately from `--max-sens`, and `deriv2`
requires a finite compile-time AD width, which the harness sets to
exactly the number of active sensitivity parameters.

**Accuracy** is measured against a CVODE reference at `atol = 1e-14`,
`rtol = 1e-12`, deliberately produced by the *other* implementation, so
cppDE is never scored against itself. The reference must also stay
strictly tighter than every benchmarked tolerance: if it matched the
tightest swept cell, CVODE would be scored against itself there, record
an error of exactly zero, and vanish from the work-precision plot as
though it could not reach that accuracy. The metric is a per-state scaled
max-norm: each state is divided by the largest value it reaches along
the reference trajectory, so a species living at `1e-9` counts as much
as one at `1e5`.

---

## The problems

### 1. The PEtab benchmark collection

Hass et al. (2019), *Benchmark problems for dynamic modeling of
intracellular processes*, **Bioinformatics** 35(17):3073–3082,
<https://doi.org/10.1093/bioinformatics/btz020>, maintained at
<https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab>.

`fetch-models.R` clones it into `benchmarks/cache/` (git-ignored; the
models keep their own licence and are never vendored into this repo).

The models are read **directly from their SBML**, not transcribed by
hand: `R/sbml.R` is a small SBML reader covering exactly what the
collection uses (compartments, species, parameters, `initialAssignment`,
`assignmentRule`, `rateRule`, `functionDefinition`, reaction kinetic
laws and MathML `piecewise`) and emits the R-syntax strings cppDE
consumes. `R/petab.R` then supplies nominal parameter values, the
condition-specific overrides, and the measurement time grid.

All 35 models parse. System sizes span **3 to 1228 states** and 5 to
4852 parameters, which is the main reason to prefer this collection over
a hand-picked set: it covers three orders of magnitude in problem size
with realistic parameter values.

### 2. Classical stiff test problems

From Hairer & Wanner, *Solving Ordinary Differential Equations II*, and
the associated IVP test set, the failure modes the biological models do
not exercise:

Robertson, HIRES, OREGO (Oregonator), E5, Van der Pol at μ = 1000,
Pollution (20 species / 25 reactions), plus two problems whose size is a
free knob: a 1-D Brusselator method-of-lines discretisation and a
FitzHugh–Nagumo chain, both with banded Jacobians for the sparse path.

---

## Modelling decisions, and why

These are the points where a PEtab problem does not map one-to-one onto
"integrate this ODE". Each is applied identically to both backends, so
the comparison stays fair, and each is recorded in the problem's `notes`
and printed by the runner, nothing is silently simplified.

**Time switches become events.** Several models gate a parameter on
time via MathML `piecewise` (a dose at `t₀`, a lockdown window, a
ramp that stops at `EGF_end`). A switch that falls *outside* the
simulated window folds to a constant. A switch *inside* it is a genuine
discontinuity, and is turned into a **timed event** on an auxiliary
state with `dx/dt = 0`: the solver stops at the switch, applies the
jump, and restarts, instead of stepping across a discontinuity and
letting the error controller discover it. The switch instant is located
by bisection. Identical `piecewise` expressions share one auxiliary
state, which matters for `Chen_MSB2009`, where the same switch appears
about a hundred times.

Where a branch is itself time-dependent, 0/1 indicator states select
between the branch expressions instead.

*Consequence:* a parameter that controls a switch *time* is dropped
from the sensitivity set, because event times are emitted as numeric
constants and the saltation term such a parameter would contribute is
therefore not in the model. Reporting no gradient is better than
reporting a wrong one.

**Sensitivity parameters.** The set is PEtab's `estimate == 1`
parameters, intersected with those that actually appear in the ODE,
so observable and noise parameters (`sd_*`, `offset_*`, `scaling_*`)
drop out on their own. Species initials are held fixed unless the model
estimates them. `--max-sens` caps the count (default 32); the cap is
reported whenever it bites.

**Pre-equilibration is not performed.** Fifteen models specify a
steady-state pre-equilibration condition. Trajectories start from the
SBML/`initialAssignment` initials instead. This changes *which*
trajectory is integrated but not how hard it is, and both backends
integrate the same one.

**Conditions.** `--conditions all` runs every experimental condition of
each model, which is the realistic unit of work for a fitting run: one
optimiser step solves them all. The default is one condition per model
(the one with the richest finite time grid). The right-hand side is
compiled once and reused across conditions, except where condition
values change how a `piecewise` folds, which is detected and handled.

**Excluded.** `Smith_BMCSystBiol2013` alone switches a `piecewise` on a
*state* rather than on time; it is reported and skipped
(`--include-excluded` runs it anyway).

---

## The figures

| figure | question it answers |
|---|---|
| `01/02-work-precision` | how much time does a given accuracy cost? (the standard ODE-solver diagram; down-and-left is better) |
| `03-speedup` | per problem, how does cppDE compare to CVODE? Bars grow from parity, on a log₂ ratio axis |
| `04-scaling` | how does cost grow with the number of states? |
| `05-sens-overhead` | what does one solve-with-gradient cost, in units of plain solves? The grey line is the *M+1* cost of finite differences |
| `06-summary` | one geometric-mean number per mode, over every matched cell |
| `07-sens2-cost` | what a full Hessian costs versus M: cppDE only, a cost curve rather than a comparison |

Ratios are always combined with a **geometric mean**, which is the
correct average for relative measures. Cells are matched on
(problem, condition, mode, tolerance) before any ratio is taken.

For work-precision diagrams use the finer tolerance sweep:

```sh
Rscript benchmarks/run-benchmarks.R --tol wp --nrep 7
```

---

## Running it somewhere else

A full sweep is long and wants an idle machine, so there are two
submitters. Both parse the problems **here** and ship them inside the
transferred workspace, so the remote side needs neither this repository
nor the PEtab collection, only cppDE with the CVODE backend, a C++
compiler, and R. Both split the work with the same cost-balanced
sharding (`R/shards.R`), so their results are comparable with each other
and with a local run.

| script | via | placement |
|---|---|---|
| `run-on-cluster.R` | `dMod2::distributedComputing()` | a SLURM array; the scheduler reserves cores and memory |
| `run-with-runbg.R` | `dMod2::runbg()` | plain ssh; one background `R CMD BATCH` per shard, placed by you |

```sh
Rscript benchmarks/run-with-runbg.R --dry-run --tier full --shards 4
Rscript benchmarks/run-with-runbg.R --submit  --machines you@box --tier full --shards 4
Rscript benchmarks/run-with-runbg.R --check   --jobname cppde_bg
Rscript benchmarks/run-with-runbg.R --log     --jobname cppde_bg
Rscript benchmarks/run-with-runbg.R --collect --jobname cppde_bg
Rscript benchmarks/run-with-runbg.R --purge   --jobname cppde_bg
```

`--collect` replays the submission's options from
`results/.jobs/<jobname>.runbg.rds`, so it needs nothing but the job
name. It merges the shards into a normal results folder (CSV, README
and figures) tagged `_runbg_`.

Without a scheduler, **placement is entirely yours**. `--machines` is a
pool of hosts and `--shards n` deals *n* jobs round-robin over it, so one
host named once with `--shards 4` gets four concurrent R processes;
`--submit` prints the resulting `jobs × workers` per host. Keep that
product inside the machine, remembering that a worker can peak in the
gigabytes while it compiles. Nothing reserves the host either, so
co-tenants land in the timings: ratios stay valid, both solvers of a
matched cell run in the same process, while absolute milliseconds are
only comparable within one host. Every row carries `node` and `cpu`.

`--submit` first probes each host over ssh for R, cppDE, the CVODE
backend, KLU and a C++ compiler, and refuses to send anything if one is
missing. runbg starts `R CMD BATCH --vanilla` through a *non-login* ssh
shell, so whatever a module system or `~/.profile` sets up has to be
reachable from `~/.bashrc`. Authentication must be non-interactive:
unlike `distributedComputing()`, `runbg()` has no `ssh_passwd`, so a key
or an agent is required, including for `--machines localhost`.

---

## Memory guards

The peak of a benchmark run is the *build*, not the solve: a large model
with sensitivities is a very large translation unit and g++ answers it in
gigabytes. With one worker per problem those compilers run concurrently,
and on a host without a scheduler nothing bounds their sum. Enough of
them at once does not merely swap, the kernel goes into reclaim and
stops answering, and the job dies with the machine.

Three guards bound that, each covering what the others cannot. They apply
to every route (local, runbg, SLURM) with the same flag names:

| flag | bounds | default |
|---|---|---|
| `--compile-slots <n>` | concurrent compilers per **host** | `min(cores, 4)` |
| `--max-compile-gb <g>` | address space of one compiler | `8` |
| `--max-worker-gb <g>` | address space of one worker | off |

`--compile-slots` is a counting semaphore under `$TMPDIR`, so it is
host-local and shared: two benchmark jobs on the same box count against
the same slots, which is the case that matters, because that is the one
no single process can see. `0` disables it.

`--max-compile-gb` is a `ulimit -v` applied through a generated
`R_MAKEVARS_USER` that wraps the toolchain R reports. A model that needs
more than the cap fails to build, is logged as a skipped model, and the
run continues, instead of the host going down with it. Together the two
put a host's build footprint at about `compile-slots × max-compile-gb`,
which `--submit` prints before it sends anything.

`--max-worker-gb` caps the R worker itself and needs the `unix` package
on the hosts; without it the run says so and continues, since the
compiler cap is the one that always applies.

For the ssh route the preflight reads each host's cores, total and
available RAM and load, and **refuses to submit** unless three things
hold on every host:

* not more workers than cores, `jobs × --bench-cores ≤ nproc`
* at least `--mem-per-worker` (6 GB) of available RAM per worker
* a build ceiling of `min(--compile-slots, workers) × --max-compile-gb`
  that stays under `--mem-budget` (0.5) of what is free

The third is the one that makes "it will not take the machine down" a
checked property rather than an assumption: it compares the guards' own
maximum against the machine they are about to run on, and refuses a run
with no `--max-compile-gb` at all. Planning for only half of free memory
is deliberate, `MemAvailable` is a snapshot of a box nothing reserves,
and the co-tenant who logs in next has to fit somewhere too. `--force`
turns any of these refusals into a warning.

```
  knecht5   R 4.5.0 | cppDE yes | CVODE yes | KLU yes | C++ yes | solve yes
            12 cores, 504 GB RAM (498 GB available), load 0.14 -> 6 job(s) x 8 = 48 worker(s)
            TOO MANY WORKERS: 48 workers on 12 cores
            build ceiling 32 GB (4 compiler(s) x 8 GB) vs 249 GB budget (50% of free)
```

Under SLURM the allocation already reserves memory, so there the guards
are a second line: what they add is that an oversized model fails to
build and is skipped, rather than the task being cancelled for exceeding
its reservation.

---

## Layout

```
benchmarks/
  fetch-models.R        download / update the PEtab collection
  run-benchmarks.R      the driver -- start here
  run-on-cluster.R      submit to SLURM     (dMod2::distributedComputing)
  run-with-runbg.R      submit over ssh     (dMod2::runbg)
  validate-models.R     check the SBML translation against the
                        simulated data shipped with the collection
  R/
    sbml.R              SBML + MathML  ->  R expression strings
    petab.R             PEtab tables   ->  benchmark problems
    problems-classic.R  the classical stiff test problems
    shards.R            problem selection + cost-balanced splitting,
                        shared by both submitters
    harness.R           compile / solve / time / score
    plots.R             the ggplot2 figures
  cache/                fetched models        (git-ignored)
  results/              CSV + figures         (git-ignored)
```

Adding a problem means appending a `bench_problem()` to
`R/problems-classic.R`; the harness and every figure pick it up
automatically.

---

## Requirements

* cppDE installed with the **CVODE backend enabled**, the run aborts
  early with instructions otherwise. `cppDE::install_libs("sundials")`
  builds SUNDIALS into a per-user cache without root.
* R packages `ggplot2`, `scales`, `xml2`, `yaml`.
* `git`, for `fetch-models.R`.

Runs compile a lot of C++. A full sweep is dominated by compilation, not
by integration; compiled models are cached per session and reused across
tolerances and conditions, and released before the next model so the run
stays inside R's limit on simultaneously loaded shared libraries.

**Run it on an idle machine.** Timings are wall-clock. Each cell is
repeated until it has run for at least 50 ms and the median of `--nrep`
such batches is reported, which removes clock granularity but not a
competing compile job. Anything else running on the box will show up in
the numbers.

### `--cores`: what parallelism is safe here

`--cores n` runs *n problems* concurrently. The granularity is
deliberate: one worker executes **all** of a problem's solvers, modes and
tolerances in sequence, so both sides of every matched cell run under the
same machine load.

That matters because a benchmark under load is not simply "the same
numbers, later". Cores on a socket share last-level cache and memory
bandwidth, SMT siblings share execution units, and, usually the largest
effect, CPUs clock down as more cores go busy, so a core that runs at
its turbo ceiling alone runs measurably slower with twenty neighbours.

Those effects are largely *uniform*, so a **ratio** between two solvers
on the same cell cancels much of them, which is why cells are never
split across workers. But "much of" is not "all of", and the suite was
measured to find out how much. Running the `tiny` tier at 1 and at 8
workers on the same machine:

| | serial | 8 workers | shift |
|---|---:|---:|---:|
| ratio, with sensitivities | 1.185 | 1.165 | −1.7 % |
| ratio, without sensitivities | 0.791 | 0.890 | **+12.5 %** |
| worst single cell, no sens | | | **96 %** |

The split is explained by duration, not by the mode: a `nosens` solve on
these problems takes 0.4–2.3 ms, and at that scale scheduler jitter under
load dominates. The `sens1` solves are one to two orders of magnitude
longer and come out stable. `bench_time()` therefore averages over
250 ms batches instead of 50 ms when `--cores > 1`, which is what buys
the short cells back.

Practical rule:

* **`--cores 1`** for any number you intend to quote, and for judging a
  change smaller than ~10 %;
* **`--cores > 1`** for a fast "did I break something" pass, treat
  differences below ~10 % as noise;
* **absolute milliseconds** are inflated either way at `--cores > 1`;
  never compare them across runs made with different `--cores` values.

The value is written to the `cores` column of `results.csv` and to
`run-info.txt`, so the two kinds of run cannot be mixed up after the
fact. Workers get `OMP_NUM_THREADS=1` so a threaded BLAS cannot
oversubscribe on top of the fan-out.

Translating the SBML models is parallelised too, unconditionally, it
happens before anything is timed.

---

## Two cppDE issues this suite surfaced

Both are worked around in `R/sbml.R` so the benchmark runs, but they are
package-level bugs, not benchmark-level ones:

1. **Parameter names that are C++ reserved words generate uncompilable
   code.** SymPy's C++ printer renames a reserved word by appending an
   underscore (`default` → `default_`), while cppDE declares the
   parameter under its original name. The generated file then fails with
   `'default_' was not declared in this scope`. This is not exotic: the
   SBML default compartment is called `default`, so it affects a large
   share of the collection. *Workaround:* `san_id()` renames such
   identifiers on ingest.

2. **A CSE temporary holding a constant subexpression is typed as the AD
   scalar.** An expression like `log(2)/tau` makes CSE emit
   `const AD _cse_t2 = cppde::log(2.0);`, and there is no
   `cppde::log(double)` overload, so the model compiles with
   `deriv = FALSE` and fails only with `deriv = TRUE`. *Workaround:*
   `fold_constants()` evaluates every symbol-free subexpression before
   codegen.

## Numerical observations worth following up

Cross-checking the two backends on the identical right-hand side at
`atol = 1e-10`, `rtol = 1e-8` puts nearly every model at a relative
agreement of `1e-8` or better. Three do not, and they are worth a look:

* `Crauste_CellSystems2017`, backends differ by `4e-3` in the states.
  The same model is also the one that misses the `simulatedData` check
  by `1.5e-3`, so the difficulty is in the integration, not the
  translation.
* `Isensee_JCB2018`, differs by O(1) *without* sensitivities while
  agreeing to `6e-10` *with* them. That direction is expected, not
  paradoxical: with sensitivities the error test also covers the
  sensitivity variables, so at the same tolerance the controller is
  forced onto smaller steps and both backends land much closer to the
  true trajectory. The finding is therefore about the plain run, at
  `rtol = 1e-8` the requested tolerance is not enough to pin this model
  down, and the two step-size controllers diverge.
* `Elowitz_Nature2000`, states agree to `5e-7`, sensitivities do not.
  The model estimates a Hill exponent `n_Hill`, and `d(x^n)/dn` involves
  `log(x)` while some states approach zero, so the sensitivity is
  genuinely ill-posed there rather than wrongly computed.
