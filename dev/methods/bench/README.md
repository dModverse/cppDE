# Benchmark run `20260901-182706_full_c1_nosens-sens1`

cppDE against SUNDIALS CVODE(S). 852 rows, 4 failed cells.

## Result

| mode | solver | time vs CVODE | rhs-evals vs CVODE | cells | problems |
|---|---|---:|---:|---:|---:|
| without sensitivities | cppDE_ndf | **1.20×** | 0.95× | 123 | 41 |
| with 1st-order sensitivities | cppDE_ndf | **1.84×** | 0.08× | 121 | 41 |

## Dense vs sparse LU

15 model(s) run a second and third time with the linear solver pinned, next to the auto-detected head-to-head.
Gain > 1 means the sparse path was faster; `chose` is what
auto-detection picked, so a gain > 1 next to `dense` is a
mis-detection.

| problem | states | backend | mode | chose | gain |
|---|---:|---|---|---|---:|
| Elowitz_Nature2000 | 8 | cppde | without sensitivities | sparse | **1.13×** |
| Elowitz_Nature2000 | 8 | cppde | with 1st-order sensitivities | sparse | **1.26×** |
| Elowitz_Nature2000 | 8 | cvode | without sensitivities | sparse | **0.94×** |
| Elowitz_Nature2000 | 8 | cvode | with 1st-order sensitivities | sparse | **1.01×** |
| Oliveira_NatCommun2021 | 10 | cppde | without sensitivities | sparse | **1.06×** |
| Oliveira_NatCommun2021 | 10 | cppde | with 1st-order sensitivities | sparse | **1.24×** |
| Oliveira_NatCommun2021 | 10 | cvode | without sensitivities | sparse | **1.09×** |
| Oliveira_NatCommun2021 | 10 | cvode | with 1st-order sensitivities | sparse | **1.01×** |
| Raia_CancerResearch2011 | 14 | cppde | without sensitivities | sparse | **1.18×** |
| Raia_CancerResearch2011 | 14 | cppde | with 1st-order sensitivities | sparse | **1.10×** |
| Raia_CancerResearch2011 | 14 | cvode | without sensitivities | sparse | **1.11×** |
| Raia_CancerResearch2011 | 14 | cvode | with 1st-order sensitivities | sparse | **1.17×** |
| Pollution | 20 | cppde | without sensitivities | sparse | **1.24×** |
| Pollution | 20 | cppde | with 1st-order sensitivities | sparse | **1.09×** |
| Pollution | 20 | cvode | without sensitivities | sparse | **1.32×** |
| Pollution | 20 | cvode | with 1st-order sensitivities | sparse | **1.01×** |
| Raimundez_PCB2020 | 22 | cppde | without sensitivities | sparse | **1.27×** |
| Raimundez_PCB2020 | 22 | cppde | with 1st-order sensitivities | sparse | **1.08×** |
| Raimundez_PCB2020 | 22 | cvode | without sensitivities | sparse | **1.14×** |
| Raimundez_PCB2020 | 22 | cvode | with 1st-order sensitivities | sparse | **1.03×** |
| Bachmann_MSB2011 | 25 | cppde | without sensitivities | sparse | **1.45×** |
| Bachmann_MSB2011 | 25 | cppde | with 1st-order sensitivities | sparse | **1.22×** |
| Bachmann_MSB2011 | 25 | cvode | without sensitivities | sparse | **1.61×** |
| Bachmann_MSB2011 | 25 | cvode | with 1st-order sensitivities | sparse | **1.03×** |
| Isensee_JCB2018 | 25 | cppde | without sensitivities | sparse | **1.33×** |
| Isensee_JCB2018 | 25 | cppde | with 1st-order sensitivities | sparse | **1.27×** |
| Isensee_JCB2018 | 25 | cvode | without sensitivities | sparse | **1.42×** |
| Isensee_JCB2018 | 25 | cvode | with 1st-order sensitivities | sparse | **1.06×** |
| Lucarelli_CellSystems2018 | 33 | cppde | without sensitivities | sparse | **1.49×** |
| Lucarelli_CellSystems2018 | 33 | cppde | with 1st-order sensitivities | sparse | **1.15×** |
| Lucarelli_CellSystems2018 | 33 | cvode | without sensitivities | sparse | **1.60×** |
| Lucarelli_CellSystems2018 | 33 | cvode | with 1st-order sensitivities | sparse | **1.05×** |
| Laske_PLOSComputBiol2019 | 34 | cppde | without sensitivities | sparse | **1.54×** |
| Laske_PLOSComputBiol2019 | 34 | cppde | with 1st-order sensitivities | sparse | **1.41×** |
| Laske_PLOSComputBiol2019 | 34 | cvode | without sensitivities | sparse | **1.68×** |
| Laske_PLOSComputBiol2019 | 34 | cvode | with 1st-order sensitivities | sparse | **1.11×** |
| Alkan_SciSignal2018 | 36 | cppde | without sensitivities | sparse | **1.29×** |
| Alkan_SciSignal2018 | 36 | cppde | with 1st-order sensitivities | sparse | **1.27×** |
| Alkan_SciSignal2018 | 36 | cvode | without sensitivities | sparse | **1.22×** |
| Alkan_SciSignal2018 | 36 | cvode | with 1st-order sensitivities | sparse | **1.00×** |
| Brusselator1D_N24 | 48 | cppde | without sensitivities | sparse | **1.57×** |
| Brusselator1D_N24 | 48 | cppde | with 1st-order sensitivities | sparse | **1.84×** |
| Brusselator1D_N24 | 48 | cvode | without sensitivities | sparse | **1.68×** |
| Brusselator1D_N24 | 48 | cvode | with 1st-order sensitivities | sparse | **1.23×** |
| FitzHughNagumo_N24 | 48 | cppde | without sensitivities | sparse | **1.89×** |
| FitzHughNagumo_N24 | 48 | cppde | with 1st-order sensitivities | sparse | **1.63×** |
| FitzHughNagumo_N24 | 48 | cvode | without sensitivities | sparse | **1.95×** |
| FitzHughNagumo_N24 | 48 | cvode | with 1st-order sensitivities | sparse | **1.21×** |
| Giordano_Nature2020 | 51 | cppde | without sensitivities | sparse | **2.15×** |
| Giordano_Nature2020 | 51 | cppde | with 1st-order sensitivities | sparse | **1.84×** |
| Giordano_Nature2020 | 51 | cvode | without sensitivities | sparse | **2.50×** |
| Giordano_Nature2020 | 51 | cvode | with 1st-order sensitivities | sparse | **1.19×** |
| Brusselator1D_N64 | 128 | cppde | without sensitivities | sparse | **2.44×** |
| Brusselator1D_N64 | 128 | cppde | with 1st-order sensitivities | sparse | **2.41×** |
| Brusselator1D_N64 | 128 | cvode | without sensitivities | sparse | **3.23×** |
| Brusselator1D_N64 | 128 | cvode | with 1st-order sensitivities | sparse | **1.73×** |
| FitzHughNagumo_N64 | 128 | cppde | without sensitivities | sparse | **3.66×** |
| FitzHughNagumo_N64 | 128 | cppde | with 1st-order sensitivities | sparse | **2.54×** |
| FitzHughNagumo_N64 | 128 | cvode | without sensitivities | sparse | **3.80×** |
| FitzHughNagumo_N64 | 128 | cvode | with 1st-order sensitivities | sparse | **1.68×** |

## Configuration

| | |
|---|---|
| tier | full |
| modes | nosens,sens1 |
| cores | 1 (serial) |
| tolerances | atol 1e-06 / rtol 1e-04; atol 1e-09 / rtol 1e-07; atol 1e-12 / rtol 1e-10 |
| repetitions | 5 |
| elapsed | 2661.7 s |
| date | 2026-09-01 19:11:35 |
| R / cppDE | 4.6.1 / 0.9.2 |
| CPU | Intel(R) Core(TM) Ultra 9 185H |
| KLU | TRUE |

Full options: `tier=full suite=all models= skip=Lang_PLOSComputBiol2024 conditions=1 modes=nosens,sens1 max-sens2=10 tol=default nrep=5 min-time= cores=1 compile-slots=1 max-compile-gb=8 max-worker-gb= shard= max-states=400 max-sens=32 min-points=25 outdir=/home/simon/Documents/Projects/dModverse/cppDE/benchmarks/results petab-root=/home/simon/Documents/Projects/dModverse/cppDE/benchmarks/cache/petab/Benchmark-Models builddir= extra-solvers=FALSE sparse-sweep=TRUE max-density=0.25 min-sweep-states=8 include-excluded=FALSE plots=TRUE quick=FALSE help=FALSE`

## Problems

| problem | source | states | params | sens | conditions | out |
|---|---|---:|---:|---:|---:|---:|
| VanDerPol_mu1000 | classic | 2 | 1 | 1 | 1 | 500 |
| Bertozzi_PNAS2020 | petab | 3 | 5 | 3 | 1 | 33 |
| Borghans_BiophysChem1997 | petab | 3 | 22 | 20 | 1 | 112 |
| OREGO | classic | 3 | 3 | 3 | 1 | 400 |
| Robertson | classic | 3 | 3 | 3 | 1 | 201 |
| Armistead_CellDeathDis2024 | petab | 4 | 12 | 10 | 1 | 25 |
| E5 | classic | 4 | 4 | 3 | 1 | 201 |
| Perelson_Science1996 | petab | 4 | 6 | 2 | 1 | 39 |
| Zhao_QuantBiol2020 | petab | 4 | 33 | 21 | 1 | 37 |
| Beer_MolBioSystems2014 | petab | 5 | 8 | 6 | 1 | 714 |
| Crauste_CellSystems2017 | petab | 5 | 13 | 12 | 1 | 31 |
| Fiedler_BMCSystBiol2016 | petab | 6 | 18 | 12 | 1 | 30 |
| Sneyd_PNAS2002 | petab | 6 | 17 | 14 | 1 | 37 |
| Bruno_JExpBot2016 | petab | 7 | 18 | 7 | 1 | 26 |
| Liu_IFACPapersOnLine2025 | petab | 7 | 8 | 7 | 1 | 25 |
| Rahman_MBS2016 | petab | 7 | 19 | 9 | 1 | 45 |
| Boehm_JProteomeRes2014 | petab | 8 | 9 | 6 | 1 | 28 |
| Elowitz_Nature2000 | petab | 8 | 19 | 18 | 1 | 59 |
| HIRES | classic | 8 | 2 | 2 | 1 | 200 |
| Okuonghae_ChaosSolitonsFractals2020 | petab | 8 | 17 | 14 | 1 | 46 |
| Weber_BMC2015 | petab | 8 | 34 | 26 | 1 | 33 |
| Brannmark_JBC2010 | petab | 9 | 19 | 14 | 1 | 32 |
| Fujita_SciSignal2010 | petab | 9 | 20 | 16 | 1 | 27 |
| Oliveira_NatCommun2021 | petab | 10 | 24 | 10 | 1 | 60 |
| Schwen_PONE2014 | petab | 11 | 15 | 13 | 1 | 32 |
| Raia_CancerResearch2011 | petab | 14 | 21 | 18 | 1 | 26 |
| Zheng_PNAS2012 | petab | 15 | 47 | 45 | 1 | 27 |
| Blasi_CellSystems2016 | petab | 16 | 10 | 8 | 1 | 25 |
| Pollution | classic | 20 | 25 | 25 | 1 | 200 |
| Raimundez_PCB2020 | petab | 22 | 79 | 57 | 1 | 28 |
| Bachmann_MSB2011 | petab | 25 | 39 | 27 | 1 | 29 |
| Isensee_JCB2018 | petab | 25 | 59 | 32 | 1 | 27 |
| Lucarelli_CellSystems2018 | petab | 33 | 107 | 64 | 1 | 30 |
| Laske_PLOSComputBiol2019 | petab | 34 | 69 | 6 | 1 | 33 |
| Alkan_SciSignal2018 | petab | 36 | 53 | 34 | 1 | 28 |
| Brusselator1D_N24 | classic | 48 | 3 | 3 | 1 | 100 |
| FitzHughNagumo_N24 | classic | 48 | 5 | 5 | 1 | 100 |
| Giordano_Nature2020 | petab | 51 | 59 | 43 | 1 | 46 |
| SalazarCavazos_MBoC2020 | petab | 75 | 27 | 6 | 1 | 27 |
| Brusselator1D_N64 | classic | 128 | 3 | 3 | 1 | 100 |
| FitzHughNagumo_N64 | classic | 128 | 5 | 5 | 1 | 100 |

## Trait coverage

A tier is only a usable regression basis while every trait still has
a representative.

| trait | problems |
|---|---:|
| stiff-extreme | 2 |
| stiff-moderate | 8 |
| oscillatory | 2 |
| relaxation | 3 |
| sparse | 6 |
| large | 3 |
| events | 4 |
| many-sens | 11 |
| few-sens | 11 |
| rational | 14 |
| transcendental | 4 |
| log-horizon | 2 |
| systems-biology | 32 |

## Skipped

- `Chen_MSB2009`: ~500 species > --max-states 400
- `Froehlich_CellSystems2018`: ~1396 species > --max-states 400
- `Smith_BMCSystBiol2013`: piecewise switches on a state variable -- EXCLUDED

## Figures

| file | shows |
|---|---|
| `01`/`02-work-precision` | achieved error against cost, one panel per problem; down-and-left is better |
| `03-speedup` | per-problem ratio against CVODE, bars growing from parity |
| `04-scaling` | cost against number of states, log-log |
| `05-sens-overhead` | gradient cost in units of plain solves; grey line is the finite-difference cost |
| `06-summary` | one geometric-mean number per mode |
| `07-sens2-cost` | Hessian cost against M (cppDE only) |
| `08-sparse-gain` | dense/sparse ratio per model and backend, bars growing from parity |
| `09-sparse-crossover` | the same ratios against system size |

---

Raw numbers are in `results.csv`: one row per problem / condition /
solver / mode / tolerance.
