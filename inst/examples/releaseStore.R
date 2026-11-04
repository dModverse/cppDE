\donttest{
if (codegenAvailable()) {
  rev <- cppODE(c(A = "-k1*A", B = "k1*A - k2*B"), derivMode = "reverse",
                modelname = "releaseStore_example")
  parms <- c(A = 1, B = 0, k1 = 0.5, k2 = 0.2)
  first <- solveODE(rev, 0:10, parms, keepStore = TRUE)
  onB <- cbind(A = 0, B = rep(1, 11))
  solveODE(rev, 0:10, parms, cotangent = onB, store = first$store)$cotangent
  ## releaseStore() frees the checkpoints; a solve given the store integrates again
  releaseStore(first$store)
  solveODE(rev, 0:10, parms, cotangent = onB, store = first$store)$cotangent
}
}
