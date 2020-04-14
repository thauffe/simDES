gen_sim_df_cor <- function(TimeSim,
                           Origin,
                           Species = 1,
                           Qtimes = NULL,
                           CovBinned = NULL,
                           Ncat = NULL,
                           StateObserved = FALSE)
{
  # TimeSim Time trajectory (vector)
  # Origin Area of origin (character)
  # CovBinned Binned covariate (vector)
  # if (any(Origin == 0))
  # {
  #   Start <- sample(2:4, 1)
  # } else
  # {
    Start <- sample(min(Origin):max(Origin), 1) + 1
  # }

  if (is.null(CovBinned))
  {
    CovBinned <- NA
  }

  SimDf <- data.frame(subject = Species, time = TimeSim, state = Start, cov = CovBinned,
                      DivA = 0, DivB = 0, DivAB = 0, Strata = 1,
                      stateSampling = NA_integer_,
                      GammaCat = ifelse(is.null(Ncat), 1, sample(1:Ncat, 1)),
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
