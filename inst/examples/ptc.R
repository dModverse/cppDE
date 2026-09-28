\donttest{
if (codegenAvailable()) {
  ## Enzyme E converts S to P, P returns to S. E + ES and S + ES + P are
  ## conserved, so the steady state is fixed by the two totals.
  f <- c(E  = "-kon * E * S + (koff + kcat) * ES",
         S  = "-kon * E * S + koff * ES + kback * P",
         ES = " kon * E * S - (koff + kcat) * ES",
         P  = " kcat * ES - kback * P")
  enz <- cppFUN(f, variables = names(f), parameters = c("kon", "koff", "kcat", "kback"),
                modelname = "ptc_enzyme", outdir = tempdir(), compile = TRUE)
  C <- rbind(Etot = c(E = 1, S = 0, ES = 1, P = 0),
             Stot = c(E = 0, S = 1, ES = 1, P = 1))
  ss <- ptc(enz, x = c(E = 1, S = 1, ES = 1, P = 1),
            parms = c(kon = 1, koff = 0.5, kcat = 0.3, kback = 0.1),
            C = C, total = c(1, 10))
  ss$x
  enz$func(E = ss$x[["E"]], S = ss$x[["S"]], ES = ss$x[["ES"]], P = ss$x[["P"]],
           kon = 1, koff = 0.5, kcat = 0.3, kback = 0.1)
}
}
