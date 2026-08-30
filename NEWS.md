# cppDE 0.7.8

* Event engine with saltation corrections, root finding and PCHIP forcings.

# cppDE 0.7.7

* Rosenbrock4 and Tsit5 single-step steppers with embedded error estimate
  and dense output.

# cppDE 0.7.6

* Variable-order BDF/NDF and Adams multistep steppers in Nordsieck form,
  with a Newton corrector.

# cppDE 0.7.5

* AD aware dense LU through LAPACK and sparse LU through KLU.

# cppDE 0.7.4

* Thread-local bump arena and contiguous tangent slab as tangent storage.

# cppDE 0.7.3

* `dual2nd<T, N>` for second order sensitivities.

# cppDE 0.7.2

* Expression templates collapse right-hand side temporaries into one fused
  chain-rule loop.

# cppDE 0.7.1

* Forward-mode dual numbers `dual<T, N>` for first order sensitivities,
  with a static and a heap allocated specialisation.

# cppDE 0.7.0

* Package skeleton and licence.
