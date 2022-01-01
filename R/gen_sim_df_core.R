gen_sim_df_cor <- function(TimeSim,
                           Origin,
                           Species = 1,
                           Qtimes = NULL,
                           CovDisBinned = NULL,
                           CovExtBinned = NULL,
                           Ncat = NULL,
                           StateObserved = FALSE)
{
  Start <- resample(min(Origin):max(Origin), 1) + 1
  if (is.null(CovDisBinned))
  {
    CovDisBinned <- NA
  }
  if (is.null(CovExtBinned))
  {
    CovExtBinned <- NA
  }
  if (is.null(Ncat))
  {
    GammaCat <- 1
  }
  else if (is.infinite(Ncat))
  {
    GammaCat <- Inf
  }
  else
  {
    GammaCat <- sample(1:Ncat, 1)
  }
  SimDf <- data.frame(subject = Species, time = TimeSim, state = Start,
                      CovDisBinned, CovExtBinned,
                      DivA = 0, DivB = 0, DivAB = 0, Strata = 1,
                      stateSampling = NA_integer_,
                      GammaCat = GammaCat,
                      rate_d12 = NA_real_, rate_d21 = NA_real_,
                      rate_e1 = NA_real_, rate_e2 = NA_real_,
                      num_d12 = 0, num_d21 = 0,
                      StateObserved = StateObserved)

  if (!is.null(Qtimes))
  {
    SimDf$Strata <- findInterval(TimeSim, Qtimes) + 1
  }
  return(SimDf)
}
