sim_sampling <- function (SimDf, Pres, Nspecies, Step, Ncat, alpha, DataInArea, TraitS, VarTraitS, CatTraitS)
{
  Idx <- ifelse(SimDf$Strata == 1, 2, 2 * SimDf$Strata)
  PresA <- Pres[Idx - 1]
  PresB <- Pres[Idx]
  if(!is.null(Ncat)) {
    if(is.infinite(Ncat)) {
      GammaRate <- rgamma(Nspecies, alpha, alpha)
      GammaRate <- GammaRate[SimDf$subject]
    }
    else {
      GammaRate <- get_gamma_rates(alpha, Ncat)
      GammaRate <- GammaRate[SimDf$GammaCat]
    }
    PresA <- PresA * GammaRate
    PresB <- PresB * GammaRate
  }
  if (!is.null(TraitS)) {
    TraitMultiSamp <- exp(colSums(VarTraitS * t(TraitS[SimDf$subject, -1, drop = FALSE])))
    PresA <- PresA * TraitMultiSamp
    PresB <- PresB * TraitMultiSamp
  }
  if (!is.null(CatTraitS)) {
    TraitMultiSamp <- rowSums(CatTraitS[SimDf$subject, , drop = FALSE])
    PresA <- PresA * TraitMultiSamp
    PresB <- PresB * TraitMultiSamp
  }
  if (!is.null(DataInArea)) {
    DataInArea2 <- DataInArea
  }
  else {
    DataInArea2 <- 0
  }
  SA <- 1 - exp(-Step * PresA)
  SB <- 1 - exp(-Step * PresB)
  SimDf$stateSampling <- 0
  # Area A
  W1 <- which(SimDf$state == 2)
  R1 <- runif(length(W1), 0, 1)
  SimDf$stateSampling[W1] <- ifelse(R1 > SA[W1], 0, 1)
  if(DataInArea2 == 2) {
    SimDf$stateSampling[W1] <- 1
  }
  # Area B
  W2 <- which(SimDf$state == 3)
  R2 <- runif(length(W2), 0, 1)
  SimDf$stateSampling[W2] <- ifelse(R2 > SB[W2], 0, 2)
  if(DataInArea2 == 1) {
    SimDf$stateSampling[W2] <- 2
  }
  # Area AB
  W3 <- which(SimDf$state == 4)
  LenW3 <- length(W3)
  A <- ifelse(runif(LenW3, 0, 1) > SA[W3], 0, 1)
  B <- ifelse(runif(LenW3, 0, 1) > SB[W3], 0, 2)
  if(DataInArea2 == 2) {
    A <- rep(1, LenW3)
  }
  if(DataInArea2 == 1) {
    B <- rep(2, LenW3)
  }
  SimDf$stateSampling[W3] <- A + B
  SimDf$stateSampling <- SimDf$stateSampling + 1
  IdxPresent <- which(SimDf$time == max(SimDf$time))
  SimDf[IdxPresent, "stateSampling"] <- SimDf[IdxPresent, "state"]
  SimDf$stateSampling[SimDf$StateObserved] <- SimDf$state[SimDf$StateObserved]
  return(SimDf)
}
