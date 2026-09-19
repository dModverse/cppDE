# Benchmark run `20260919-203447_full_cluster_c12_nosens-sens1-sens2`

cppDE against SUNDIALS CVODE(S). 1083 rows, 6 failed cells.

> Run with 12 workers: **absolute times are inflated** by shared cache and turbo clocking. The ratios below are the quantity to read; treat differences under ~10 % as noise.

## Result

| mode | solver | time vs CVODE | rhs-evals vs CVODE | cells | problems |
|---|---|---:|---:|---:|---:|
| without sensitivities | cppDE_ndf | **1.15×** | 0.94× | 129 | 43 |
| with 1st-order sensitivities | cppDE_ndf | **1.36×** | 0.10× | 124 | 42 |

## Dense vs sparse LU

17 model(s) run a second and third time with the linear solver pinned, next to the auto-detected head-to-head.
Gain > 1 means the sparse path was faster; `chose` is what
auto-detection picked, so a gain > 1 next to `dense` is a
mis-detection.

| problem | states | backend | mode | chose | gain |
|---|---:|---|---|---|---:|
| Elowitz_Nature2000 | 8 | cppde | without sensitivities | sparse | **1.13×** |
| Elowitz_Nature2000 | 8 | cppde | with 1st-order sensitivities | sparse | **1.08×** |
| Elowitz_Nature2000 | 8 | cvode | without sensitivities | sparse | **1.06×** |
| Elowitz_Nature2000 | 8 | cvode | with 1st-order sensitivities | sparse | **1.08×** |
| Oliveira_NatCommun2021 | 10 | cppde | without sensitivities | sparse | **1.17×** |
| Oliveira_NatCommun2021 | 10 | cppde | with 1st-order sensitivities | sparse | **1.18×** |
| Oliveira_NatCommun2021 | 10 | cvode | without sensitivities | sparse | **1.11×** |
| Oliveira_NatCommun2021 | 10 | cvode | with 1st-order sensitivities | sparse | **1.13×** |
| Raia_CancerResearch2011 | 14 | cppde | without sensitivities | sparse | **1.24×** |
| Raia_CancerResearch2011 | 14 | cppde | with 1st-order sensitivities | sparse | **1.15×** |
| Raia_CancerResearch2011 | 14 | cvode | without sensitivities | sparse | **1.21×** |
| Raia_CancerResearch2011 | 14 | cvode | with 1st-order sensitivities | sparse | **1.22×** |
| Pollution | 20 | cppde | without sensitivities | sparse | **1.38×** |
| Pollution | 20 | cppde | with 1st-order sensitivities | sparse | **1.12×** |
| Pollution | 20 | cvode | without sensitivities | sparse | **1.35×** |
| Pollution | 20 | cvode | with 1st-order sensitivities | sparse | **1.22×** |
| Raimundez_PCB2020 | 22 | cppde | without sensitivities | sparse | **1.38×** |
| Raimundez_PCB2020 | 22 | cppde | with 1st-order sensitivities | sparse | **1.21×** |
| Raimundez_PCB2020 | 22 | cvode | without sensitivities | sparse | **1.42×** |
| Raimundez_PCB2020 | 22 | cvode | with 1st-order sensitivities | sparse | **1.23×** |
| Bachmann_MSB2011 | 25 | cppde | without sensitivities | sparse | **1.56×** |
| Bachmann_MSB2011 | 25 | cppde | with 1st-order sensitivities | sparse | **1.28×** |
| Bachmann_MSB2011 | 25 | cvode | without sensitivities | sparse | **1.51×** |
| Bachmann_MSB2011 | 25 | cvode | with 1st-order sensitivities | sparse | **1.31×** |
| Isensee_JCB2018 | 25 | cppde | without sensitivities | sparse | **1.33×** |
| Isensee_JCB2018 | 25 | cppde | with 1st-order sensitivities | sparse | **1.28×** |
| Isensee_JCB2018 | 25 | cvode | without sensitivities | sparse | **1.37×** |
| Isensee_JCB2018 | 25 | cvode | with 1st-order sensitivities | sparse | **1.18×** |
| Lucarelli_CellSystems2018 | 33 | cppde | without sensitivities | sparse | **1.68×** |
| Lucarelli_CellSystems2018 | 33 | cppde | with 1st-order sensitivities | sparse | **1.31×** |
| Lucarelli_CellSystems2018 | 33 | cvode | without sensitivities | sparse | **1.76×** |
| Lucarelli_CellSystems2018 | 33 | cvode | with 1st-order sensitivities | sparse | **1.37×** |
| Laske_PLOSComputBiol2019 | 34 | cppde | without sensitivities | sparse | **1.74×** |
| Laske_PLOSComputBiol2019 | 34 | cppde | with 1st-order sensitivities | sparse | **1.33×** |
| Laske_PLOSComputBiol2019 | 34 | cvode | without sensitivities | sparse | **1.98×** |
| Laske_PLOSComputBiol2019 | 34 | cvode | with 1st-order sensitivities | sparse | **1.42×** |
| Alkan_SciSignal2018 | 36 | cppde | without sensitivities | sparse | **1.45×** |
| Alkan_SciSignal2018 | 36 | cppde | with 1st-order sensitivities | sparse | **1.36×** |
| Alkan_SciSignal2018 | 36 | cvode | without sensitivities | sparse | **1.33×** |
| Alkan_SciSignal2018 | 36 | cvode | with 1st-order sensitivities | sparse | **1.26×** |
| Brusselator1D_N24 | 48 | cppde | without sensitivities | sparse | **1.78×** |
| Brusselator1D_N24 | 48 | cppde | with 1st-order sensitivities | sparse | **1.70×** |
| Brusselator1D_N24 | 48 | cvode | without sensitivities | sparse | **1.83×** |
| Brusselator1D_N24 | 48 | cvode | with 1st-order sensitivities | sparse | **1.48×** |
| FitzHughNagumo_N24 | 48 | cppde | without sensitivities | sparse | **1.82×** |
| FitzHughNagumo_N24 | 48 | cppde | with 1st-order sensitivities | sparse | **1.74×** |
| FitzHughNagumo_N24 | 48 | cvode | without sensitivities | sparse | **1.86×** |
| FitzHughNagumo_N24 | 48 | cvode | with 1st-order sensitivities | sparse | **1.32×** |
| Giordano_Nature2020 | 51 | cppde | without sensitivities | sparse | **2.50×** |
| Giordano_Nature2020 | 51 | cppde | with 1st-order sensitivities | sparse | **2.06×** |
| Giordano_Nature2020 | 51 | cvode | without sensitivities | sparse | **2.73×** |
| Giordano_Nature2020 | 51 | cvode | with 1st-order sensitivities | sparse | **1.62×** |
| Lang_PLOSComputBiol2024 | 124 | cppde | without sensitivities | sparse | **2.54×** |
| Lang_PLOSComputBiol2024 | 124 | cppde | with 1st-order sensitivities | sparse | **3.91×** |
| Lang_PLOSComputBiol2024 | 124 | cvode | without sensitivities | sparse | **2.69×** |
| Lang_PLOSComputBiol2024 | 124 | cvode | with 1st-order sensitivities | sparse | **2.24×** |
| Brusselator1D_N64 | 128 | cppde | without sensitivities | sparse | **3.35×** |
| Brusselator1D_N64 | 128 | cppde | with 1st-order sensitivities | sparse | **3.28×** |
| Brusselator1D_N64 | 128 | cvode | without sensitivities | sparse | **3.59×** |
| Brusselator1D_N64 | 128 | cvode | with 1st-order sensitivities | sparse | **2.34×** |
| FitzHughNagumo_N64 | 128 | cppde | without sensitivities | sparse | **3.29×** |
| FitzHughNagumo_N64 | 128 | cppde | with 1st-order sensitivities | sparse | **3.36×** |
| FitzHughNagumo_N64 | 128 | cvode | without sensitivities | sparse | **3.18×** |
| FitzHughNagumo_N64 | 128 | cvode | with 1st-order sensitivities | sparse | **1.73×** |
| Chen_MSB2009 | 504 | cppde | without sensitivities | sparse | **3.16×** |
| Chen_MSB2009 | 504 | cppde | with 1st-order sensitivities | sparse | **3.46×** |
| Chen_MSB2009 | 504 | cvode | without sensitivities | sparse | **3.47×** |

**Second order**: the Hessian of the summed outputs, cppDE only;
CVODES has no second-order sensitivities, so this is a cost, not a comparison:

| problem | mode | M | plain [ms] | Hessian [ms] | factor |
|---|---|---:|---:|---:|---:|
| Armistead_CellDeathDis2024 | forward-forward | 10 | 0.21 | 6.9 | 32× |
| Armistead_CellDeathDis2024 | forward-reverse | 10 | 0.21 | 1.9 | 9× |
| Bachmann_MSB2011 | forward-reverse | 27 | 0.50 | 31.2 | 63× |
| Beer_MolBioSystems2014 | forward-forward | 6 | 0.29 | 6.2 | 21× |
| Beer_MolBioSystems2014 | forward-reverse | 6 | 0.29 | 2.2 | 7× |
| Bertozzi_PNAS2020 | forward-forward | 3 | 0.15 | 0.3 | 2× |
| Bertozzi_PNAS2020 | forward-reverse | 3 | 0.15 | 0.3 | 2× |
| Blasi_CellSystems2016 | forward-reverse | 8 | 0.32 | 4.6 | 14× |
| Boehm_JProteomeRes2014 | forward-forward | 6 | 0.34 | 13.2 | 38× |
| Boehm_JProteomeRes2014 | forward-reverse | 6 | 0.34 | 4.0 | 12× |
| Borghans_BiophysChem1997 | forward-forward | 10 | 1.26 | 360.0 | 286× |
| Borghans_BiophysChem1997 | forward-reverse | 10 | 1.26 | 26.5 | 21× |
| Brannmark_JBC2010 | forward-forward | 10 | 0.39 | 59.6 | 152× |
| Brannmark_JBC2010 | forward-reverse | 10 | 0.39 | 9.9 | 25× |
| Bruno_JExpBot2016 | forward-forward | 7 | 0.17 | 2.3 | 14× |
| Bruno_JExpBot2016 | forward-reverse | 7 | 0.17 | 0.6 | 4× |
| Crauste_CellSystems2017 | forward-forward | 10 | 0.36 | 72.3 | 198× |
| Crauste_CellSystems2017 | forward-reverse | 10 | 0.36 | 5.7 | 16× |
| E5 | forward-forward | 3 | 1.06 | 17.8 | 17× |
| E5 | forward-reverse | 3 | 1.06 | 18.4 | 17× |
| Elowitz_Nature2000 | forward-forward | 10 | 0.78 | 200.0 | 256× |
| Elowitz_Nature2000 | forward-reverse | 10 | 0.78 | 16.8 | 21× |
| Fiedler_BMCSystBiol2016 | forward-forward | 10 | 0.25 | 34.9 | 138× |
| Fiedler_BMCSystBiol2016 | forward-reverse | 10 | 0.25 | 3.1 | 12× |
| Fujita_SciSignal2010 | forward-forward | 10 | 0.28 | 40.8 | 145× |
| Fujita_SciSignal2010 | forward-reverse | 10 | 0.28 | 3.7 | 13× |
| HIRES | forward-forward | 2 | 0.42 | 5.3 | 13× |
| HIRES | forward-reverse | 2 | 0.42 | 8.4 | 20× |
| Isensee_JCB2018 | forward-reverse | 32 | 0.29 | 10.0 | 35× |
| Liu_IFACPapersOnLine2025 | forward-forward | 7 | 0.17 | 4.5 | 26× |
| Liu_IFACPapersOnLine2025 | forward-reverse | 7 | 0.17 | 0.8 | 5× |
| Okuonghae_ChaosSolitonsFractals2020 | forward-forward | 10 | 0.20 | 19.9 | 100× |
| Okuonghae_ChaosSolitonsFractals2020 | forward-reverse | 10 | 0.20 | 1.8 | 9× |
| Oliveira_NatCommun2021 | forward-forward | 10 | 0.28 | 8.4 | 30× |
| Oliveira_NatCommun2021 | forward-reverse | 10 | 0.28 | 2.9 | 10× |
| OREGO | forward-forward | 3 | 1.15 | 26.6 | 23× |
| OREGO | forward-reverse | 3 | 1.15 | 36.0 | 31× |
| Perelson_Science1996 | forward-forward | 2 | 0.21 | 1.3 | 6× |
| Perelson_Science1996 | forward-reverse | 2 | 0.21 | 1.2 | 6× |
| Pollution | forward-reverse | 25 | 0.40 | 16.9 | 42× |
| Rahman_MBS2016 | forward-forward | 9 | 0.23 | 34.4 | 148× |
| Rahman_MBS2016 | forward-reverse | 9 | 0.23 | 2.5 | 11× |
| Raia_CancerResearch2011 | forward-reverse | 18 | 0.51 | 37.5 | 73× |
| Raimundez_PCB2020 | forward-reverse | 32 | 0.39 | 28.0 | 71× |
| Robertson | forward-forward | 3 | 0.46 | 5.9 | 13× |
| Robertson | forward-reverse | 3 | 0.46 | 5.3 | 12× |
| Schwen_PONE2014 | forward-reverse | 13 | 0.34 | 6.6 | 19× |
| Sneyd_PNAS2002 | forward-forward | 10 | 0.26 | 68.5 | 265× |
| Sneyd_PNAS2002 | forward-reverse | 10 | 0.26 | 8.0 | 31× |
| VanDerPol_mu1000 | forward-forward | 1 | 0.99 | 7.4 | 8× |
| VanDerPol_mu1000 | forward-reverse | 1 | 0.99 | 13.4 | 14× |
| Weber_BMC2015 | forward-forward | 10 | 0.45 | 98.7 | 220× |
| Weber_BMC2015 | forward-reverse | 10 | 0.45 | 10.3 | 23× |
| Zhao_QuantBiol2020 | forward-forward | 10 | 0.22 | 18.7 | 87× |
| Zhao_QuantBiol2020 | forward-reverse | 10 | 0.22 | 1.7 | 8× |
| Zheng_PNAS2012 | forward-reverse | 32 | 0.34 | 13.0 | 38× |

## Configuration

| | |
|---|---|
| tier | full |
| modes | nosens,sens1,sens2 |
| cores | 12 (parallel) |
| repetitions | 5 |
| date | 2026-09-19 20:34:47 |

Full options: `jobname=cppde_bench_0919 machine=helix partition=cpu-single cores=64 bench-cores=12 mem-per-core=2 walltime=08:00:00 shards=8 tier=full modes=nosens,sens1,sens2 conditions=1 max-states=520 max-sens=32 max-sens2=10 max-sens2-fr=32 max-states-sens2=30 max-states-ff=10 reverse-from=120 models= tol=default nrep=5 min-time=0.25 sparse-sweep=TRUE max-density=0.25 min-sweep-states=8 compile-slots=4 max-compile-gb=8 max-worker-gb= outdir=/home/simon/Documents/Projects/dModverse/cppDE/benchmarks/results petab-root=/home/simon/Documents/Projects/dModverse/cppDE/benchmarks/cache/petab/Benchmark-Models ssh-passwd= libs=~/R/lib-reverseAD dry-run=FALSE submit=FALSE collect=TRUE help=FALSE`

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
