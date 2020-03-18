sim_sampling <- function(SimDf, Pres, Step, Ncat, alpha, DataInArea)
{
  Idx <- ifelse(SimDf$Strata == 1, 2, 2 * SimDf$Strata)
  PresA <- Pres[Idx - 1]
  PresB <- Pres[Idx]
  if(!is.null(Ncat))
  {
    GammaRate <- get_gamma_rates(alpha, Ncat)
    GammaRate <- GammaRate[SimDf$GammaCat]
    PresA <- PresA * GammaRate
    PresB <- PresB * GammaRate
  }
  SA <- exp(-Step * PresA)
  SB <- exp(-Step * PresB)
  SimDf$stateSampling <- 0
  # Area A
  W1 <- which(SimDf$state == 2)
  R1 <- runif(length(W1), 0, 1)
  SimDf$stateSampling[W1] <- ifelse(R1 < SA[W1], 0, 1)
  # Area B
  W2 <- which(SimDf$state == 3)
  R2 <- runif(length(W2), 0, 1)
  SimDf$stateSampling[W2] <- ifelse(R2 < SB[W2], 0, 2)
  # Area AB
  W3 <- which(SimDf$state == 4)
  LenW3 <- length(W3)
  A <- ifelse(runif(LenW3, 0, 1) < SA[W3], 0, 1)
  B <- ifelse(runif(LenW3, 0, 1) < SB[W3], 0, 2)
  SimDf$stateSampling[W3] <- A + B
  SimDf$stateSampling <- SimDf$stateSampling + 1
  IdxPresent <- which(SimDf$time == max(SimDf$time))
  SimDf[IdxPresent, "stateSampling"] <- SimDf[IdxPresent, "state"]
  return(SimDf)
}
