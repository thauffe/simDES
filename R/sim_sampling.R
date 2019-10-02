sim_sampling <- function(SimDf, Pres, Step, Ncat, alpha)
{
  Idx <- ifelse(SimDf$Strata == 1, 2, 2 * SimDf$Strata)
  PresA <- Pres[Idx - 1]
  PresB <- Pres[Idx]
  PresList <- list()
  if(is.null(Ncat))
  {
    MaxCat <- 1
    PresA <- matrix(PresA, nrow = 1)
    PresB <- matrix(PresB, nrow = 1)
  } else
  {
    PresA <- sapply(PresA, function(x) x * get_gamma_rates(alpha, Ncat))
    PresB <- sapply(PresB, function(x) x * get_gamma_rates(alpha, Ncat))
    MaxCat <- Ncat
  }
  SA <- exp(-Step * PresA)
  SB <- exp(-Step * PresB)
  PresProb <- sapply(1:ncol(SA), function(x)
    get_pres_prob(SA[SimDf$GammaCat[x], x], SB[SimDf$GammaCat[x], x])[SimDf$state[x], ],
    simplify = FALSE)
  SimDf$stateSampling <- lapply(PresProb, function(x) sample(1:4, 1, prob = x))
  IdxPresent <- which(SimDf$time == max(SimDf$time))
  SimDf[IdxPresent, "stateSampling"] <- SimDf[IdxPresent, "state"]
  return(SimDf)
}
