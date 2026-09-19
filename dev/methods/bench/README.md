# Benchmark run `20260919-164334_full_cluster_c12_nosens-sens1-sens2`

cppDE against SUNDIALS CVODE(S). 1083 rows, 15 failed cells.

> Run with 12 workers: **absolute times are inflated** by shared cache and turbo clocking. The ratios below are the quantity to read; treat differences under ~10 % as noise.

## Result

| mode | solver | time vs CVODE | rhs-evals vs CVODE | cells | problems |
|---|---|---:|---:|---:|---:|
| without sensitivities | cppDE_ndf | **1.17×** | 0.94× | 129 | 43 |
| with 1st-order sensitivities | cppDE_ndf | **1.28×** | 0.09× | 121 | 41 |

## Dense vs sparse LU

17 model(s) run a second and third time with the linear solver pinned, next to the auto-detected head-to-head.
Gain > 1 means the sparse path was faster; `chose` is what
auto-detection picked, so a gain > 1 next to `dense` is a
mis-detection.

| problem | states | backend | mode | chose | gain |
|---|---:|---|---|---|---:|
| Elowitz_Nature2000 | 8 | cppde | without sensitivities | sparse | **1.09×** |
| Elowitz_Nature2000 | 8 | cppde | with 1st-order sensitivities | sparse | **1.07×** |
| Elowitz_Nature2000 | 8 | cvode | without sensitivities | sparse | **1.09×** |
| Elowitz_Nature2000 | 8 | cvode | with 1st-order sensitivities | sparse | **1.08×** |
| Oliveira_NatCommun2021 | 10 | cppde | without sensitivities | sparse | **1.11×** |
| Oliveira_NatCommun2021 | 10 | cppde | with 1st-order sensitivities | sparse | **1.18×** |
| Oliveira_NatCommun2021 | 10 | cvode | without sensitivities | sparse | **1.08×** |
| Oliveira_NatCommun2021 | 10 | cvode | with 1st-order sensitivities | sparse | **1.12×** |
| Raia_CancerResearch2011 | 14 | cppde | without sensitivities | sparse | **1.19×** |
| Raia_CancerResearch2011 | 14 | cppde | with 1st-order sensitivities | sparse | **1.15×** |
| Raia_CancerResearch2011 | 14 | cvode | without sensitivities | sparse | **1.18×** |
| Raia_CancerResearch2011 | 14 | cvode | with 1st-order sensitivities | sparse | **1.22×** |
| Pollution | 20 | cppde | without sensitivities | sparse | **1.22×** |
| Pollution | 20 | cppde | with 1st-order sensitivities | sparse | **1.15×** |
| Pollution | 20 | cvode | without sensitivities | sparse | **1.18×** |
| Pollution | 20 | cvode | with 1st-order sensitivities | sparse | **1.23×** |
| Raimundez_PCB2020 | 22 | cppde | without sensitivities | sparse | **1.30×** |
| Raimundez_PCB2020 | 22 | cppde | with 1st-order sensitivities | sparse | **1.17×** |
| Raimundez_PCB2020 | 22 | cvode | without sensitivities | sparse | **1.32×** |
| Raimundez_PCB2020 | 22 | cvode | with 1st-order sensitivities | sparse | **1.20×** |
| Bachmann_MSB2011 | 25 | cppde | without sensitivities | sparse | **1.60×** |
| Bachmann_MSB2011 | 25 | cppde | with 1st-order sensitivities | sparse | **1.26×** |
| Bachmann_MSB2011 | 25 | cvode | without sensitivities | sparse | **1.46×** |
| Bachmann_MSB2011 | 25 | cvode | with 1st-order sensitivities | sparse | **1.30×** |
| Isensee_JCB2018 | 25 | cppde | without sensitivities | sparse | **1.32×** |
| Isensee_JCB2018 | 25 | cppde | with 1st-order sensitivities | sparse | **1.26×** |
| Isensee_JCB2018 | 25 | cvode | without sensitivities | sparse | **1.35×** |
| Isensee_JCB2018 | 25 | cvode | with 1st-order sensitivities | sparse | **1.19×** |
| Lucarelli_CellSystems2018 | 33 | cppde | without sensitivities | sparse | **1.57×** |
| Lucarelli_CellSystems2018 | 33 | cppde | with 1st-order sensitivities | sparse | **1.29×** |
| Lucarelli_CellSystems2018 | 33 | cvode | without sensitivities | sparse | **1.75×** |
| Lucarelli_CellSystems2018 | 33 | cvode | with 1st-order sensitivities | sparse | **1.38×** |
| Laske_PLOSComputBiol2019 | 34 | cppde | without sensitivities | sparse | **1.74×** |
| Laske_PLOSComputBiol2019 | 34 | cppde | with 1st-order sensitivities | sparse | **1.36×** |
| Laske_PLOSComputBiol2019 | 34 | cvode | without sensitivities | sparse | **1.98×** |
| Laske_PLOSComputBiol2019 | 34 | cvode | with 1st-order sensitivities | sparse | **1.45×** |
| Alkan_SciSignal2018 | 36 | cppde | without sensitivities | sparse | **1.31×** |
| Alkan_SciSignal2018 | 36 | cppde | with 1st-order sensitivities | sparse | **1.36×** |
| Alkan_SciSignal2018 | 36 | cvode | without sensitivities | sparse | **1.36×** |
| Alkan_SciSignal2018 | 36 | cvode | with 1st-order sensitivities | sparse | **1.27×** |
| Brusselator1D_N24 | 48 | cppde | without sensitivities | sparse | **1.70×** |
| Brusselator1D_N24 | 48 | cppde | with 1st-order sensitivities | sparse | **1.67×** |
| Brusselator1D_N24 | 48 | cvode | without sensitivities | sparse | **1.80×** |
| Brusselator1D_N24 | 48 | cvode | with 1st-order sensitivities | sparse | **1.49×** |
| FitzHughNagumo_N24 | 48 | cppde | without sensitivities | sparse | **1.84×** |
| FitzHughNagumo_N24 | 48 | cppde | with 1st-order sensitivities | sparse | **1.73×** |
| FitzHughNagumo_N24 | 48 | cvode | without sensitivities | sparse | **1.82×** |
| FitzHughNagumo_N24 | 48 | cvode | with 1st-order sensitivities | sparse | **1.32×** |
| Giordano_Nature2020 | 51 | cppde | without sensitivities | sparse | **2.57×** |
| Giordano_Nature2020 | 51 | cppde | with 1st-order sensitivities | sparse | **2.09×** |
| Giordano_Nature2020 | 51 | cvode | without sensitivities | sparse | **2.73×** |
| Giordano_Nature2020 | 51 | cvode | with 1st-order sensitivities | sparse | **1.60×** |
| Lang_PLOSComputBiol2024 | 124 | cppde | without sensitivities | sparse | **2.56×** |
| Lang_PLOSComputBiol2024 | 124 | cppde | with 1st-order sensitivities | sparse | **3.91×** |
| Lang_PLOSComputBiol2024 | 124 | cvode | without sensitivities | sparse | **2.69×** |
| Brusselator1D_N64 | 128 | cppde | without sensitivities | sparse | **3.28×** |
| Brusselator1D_N64 | 128 | cppde | with 1st-order sensitivities | sparse | **3.31×** |
| Brusselator1D_N64 | 128 | cvode | without sensitivities | sparse | **3.59×** |
| Brusselator1D_N64 | 128 | cvode | with 1st-order sensitivities | sparse | **2.27×** |
| FitzHughNagumo_N64 | 128 | cppde | without sensitivities | sparse | **3.17×** |
| FitzHughNagumo_N64 | 128 | cppde | with 1st-order sensitivities | sparse | **3.38×** |
| FitzHughNagumo_N64 | 128 | cvode | without sensitivities | sparse | **3.18×** |
| FitzHughNagumo_N64 | 128 | cvode | with 1st-order sensitivities | sparse | **1.74×** |
| Chen_MSB2009 | 504 | cppde | without sensitivities | sparse | **3.14×** |
| Chen_MSB2009 | 504 | cppde | with 1st-order sensitivities | sparse | **3.44×** |
| Chen_MSB2009 | 504 | cvode | without sensitivities | sparse | **3.43×** |

**Second order**: the Hessian of the summed outputs, cppDE only;
CVODES has no second-order sensitivities, so this is a cost, not a comparison:

| problem | mode | M | plain [ms] | Hessian [ms] | factor |
|---|---|---:|---:|---:|---:|
| Armistead_CellDeathDis2024 | forward-forward | 10 | 0.21 | 7.0 | 33× |
| Armistead_CellDeathDis2024 | forward-reverse | 10 | 0.21 | 1.9 | 9× |
| Bachmann_MSB2011 | forward-reverse | 27 | 0.48 | 30.8 | 64× |
| Beer_MolBioSystems2014 | forward-forward | 6 | 0.29 | 5.4 | 18× |
| Beer_MolBioSystems2014 | forward-reverse | 6 | 0.29 | 2.1 | 7× |
| Bertozzi_PNAS2020 | forward-forward | 3 | 0.15 | 0.3 | 2× |
| Bertozzi_PNAS2020 | forward-reverse | 3 | 0.15 | 0.3 | 2× |
| Blasi_CellSystems2016 | forward-reverse | 8 | 0.32 | 4.6 | 14× |
| Boehm_JProteomeRes2014 | forward-forward | 6 | 0.33 | 12.6 | 39× |
| Boehm_JProteomeRes2014 | forward-reverse | 6 | 0.33 | 4.0 | 12× |
| Borghans_BiophysChem1997 | forward-forward | 10 | 1.34 | 364.0 | 271× |
| Borghans_BiophysChem1997 | forward-reverse | 10 | 1.34 | 26.5 | 20× |
| Brannmark_JBC2010 | forward-forward | 10 | 0.40 | 59.2 | 149× |
| Brannmark_JBC2010 | forward-reverse | 10 | 0.40 | 9.9 | 25× |
| Bruno_JExpBot2016 | forward-forward | 7 | 0.16 | 2.3 | 14× |
| Bruno_JExpBot2016 | forward-reverse | 7 | 0.16 | 0.7 | 4× |
| Crauste_CellSystems2017 | forward-forward | 10 | 0.36 | 72.3 | 203× |
| Crauste_CellSystems2017 | forward-reverse | 10 | 0.36 | 5.6 | 16× |
| E5 | forward-forward | 3 | 1.00 | 18.3 | 18× |
| E5 | forward-reverse | 3 | 1.00 | 18.3 | 18× |
| Elowitz_Nature2000 | forward-forward | 10 | 0.72 | 219.5 | 303× |
| Elowitz_Nature2000 | forward-reverse | 10 | 0.72 | 16.8 | 23× |
| Fiedler_BMCSystBiol2016 | forward-forward | 10 | 0.25 | 31.8 | 128× |
| Fiedler_BMCSystBiol2016 | forward-reverse | 10 | 0.25 | 3.1 | 12× |
| Fujita_SciSignal2010 | forward-forward | 10 | 0.28 | 38.5 | 138× |
| Fujita_SciSignal2010 | forward-reverse | 10 | 0.28 | 3.8 | 14× |
| HIRES | forward-forward | 2 | 0.41 | 5.3 | 13× |
| HIRES | forward-reverse | 2 | 0.41 | 8.3 | 20× |
| Isensee_JCB2018 | forward-reverse | 32 | 0.28 | 9.9 | 35× |
| Liu_IFACPapersOnLine2025 | forward-forward | 7 | 0.17 | 4.5 | 27× |
| Liu_IFACPapersOnLine2025 | forward-reverse | 7 | 0.17 | 0.8 | 5× |
| Okuonghae_ChaosSolitonsFractals2020 | forward-forward | 10 | 0.21 | 20.2 | 95× |
| Okuonghae_ChaosSolitonsFractals2020 | forward-reverse | 10 | 0.21 | 1.7 | 8× |
| Oliveira_NatCommun2021 | forward-forward | 10 | 0.28 | 8.9 | 32× |
| Oliveira_NatCommun2021 | forward-reverse | 10 | 0.28 | 2.9 | 10× |
| OREGO | forward-forward | 3 | 1.17 | 27.1 | 23× |
| OREGO | forward-reverse | 3 | 1.17 | 36.1 | 31× |
| Perelson_Science1996 | forward-forward | 2 | 0.21 | 1.3 | 6× |
| Perelson_Science1996 | forward-reverse | 2 | 0.21 | 1.1 | 5× |
| Pollution | forward-reverse | 25 | 0.44 | 17.3 | 39× |
| Rahman_MBS2016 | forward-forward | 9 | 0.23 | 35.4 | 154× |
| Rahman_MBS2016 | forward-reverse | 9 | 0.23 | 2.5 | 11× |
| Raia_CancerResearch2011 | forward-reverse | 18 | 0.51 | 37.9 | 74× |
| Raimundez_PCB2020 | forward-reverse | 32 | 0.38 | 28.0 | 73× |
| Robertson | forward-forward | 3 | 0.45 | 5.8 | 13× |
| Robertson | forward-reverse | 3 | 0.45 | 5.3 | 12× |
| Schwen_PONE2014 | forward-reverse | 13 | 0.35 | 6.6 | 19× |
| Sneyd_PNAS2002 | forward-forward | 10 | 0.26 | 67.3 | 256× |
| Sneyd_PNAS2002 | forward-reverse | 10 | 0.26 | 7.9 | 30× |
| VanDerPol_mu1000 | forward-forward | 1 | 1.00 | 7.5 | 7× |
| VanDerPol_mu1000 | forward-reverse | 1 | 1.00 | 13.2 | 13× |
| Weber_BMC2015 | forward-forward | 10 | 0.43 | 107.0 | 247× |
| Weber_BMC2015 | forward-reverse | 10 | 0.43 | 10.4 | 24× |
| Zhao_QuantBiol2020 | forward-forward | 10 | 0.19 | 18.7 | 101× |
| Zhao_QuantBiol2020 | forward-reverse | 10 | 0.19 | 1.8 | 10× |
| Zheng_PNAS2012 | forward-reverse | 32 | 0.34 | 12.9 | 38× |

## Configuration

| | |
|---|---|
| tier | full |
| modes | nosens,sens1,sens2 |
| cores | 12 (parallel) |
| repetitions | 5 |
| date | 2026-09-19 16:43:34 |

Full options: `jobname=cppde_bench_0918 machine=helix partition=cpu-single cores=64 bench-cores=12 mem-per-core=2 walltime=08:00:00 shards=8 tier=full modes=nosens,sens1,sens2 conditions=1 max-states=520 max-sens=32 max-sens2=10 max-sens2-fr=32 max-states-sens2=30 max-states-ff=10 reverse-from=120 tol=default nrep=5 min-time=0.25 sparse-sweep=TRUE max-density=0.25 min-sweep-states=8 compile-slots=4 max-compile-gb=8 max-worker-gb= outdir=/home/simon/Documents/Projects/dModverse/cppDE/benchmarks/results petab-root=/home/simon/Documents/Projects/dModverse/cppDE/benchmarks/cache/petab/Benchmark-Models ssh-passwd= libs=~/R/lib-reverseAD dry-run=FALSE submit=FALSE collect=TRUE help=FALSE`

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
| Zheng_PNAS2012 | petab | 15 | 47 | 32 | 1 | 27 |
| Blasi_CellSystems2016 | petab | 16 | 10 | 8 | 1 | 25 |
| Pollution | classic | 20 | 25 | 25 | 1 | 200 |
| Raimundez_PCB2020 | petab | 22 | 79 | 32 | 1 | 28 |
| Bachmann_MSB2011 | petab | 25 | 39 | 27 | 1 | 29 |
| Isensee_JCB2018 | petab | 25 | 59 | 32 | 1 | 27 |
| Lucarelli_CellSystems2018 | petab | 33 | 107 | 32 | 1 | 30 |
| Laske_PLOSComputBiol2019 | petab | 34 | 69 | 6 | 1 | 33 |
| Alkan_SciSignal2018 | petab | 36 | 53 | 32 | 1 | 28 |
| Brusselator1D_N24 | classic | 48 | 3 | 3 | 1 | 100 |
| FitzHughNagumo_N24 | classic | 48 | 5 | 5 | 1 | 100 |
| Giordano_Nature2020 | petab | 51 | 59 | 32 | 1 | 46 |
| SalazarCavazos_MBoC2020 | petab | 75 | 27 | 6 | 1 | 27 |
| Lang_PLOSComputBiol2024 | petab | 124 | 218 | 164 | 1 | 600 |
| Brusselator1D_N64 | classic | 128 | 3 | 3 | 1 | 100 |
| FitzHughNagumo_N64 | classic | 128 | 5 | 5 | 1 | 100 |
| Chen_MSB2009 | petab | 504 | 189 | 152 | 1 | 32 |

## Skipped

_none_

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
