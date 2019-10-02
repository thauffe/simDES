gen_sim_df_cor <- function(TimeSim,
                           Origin,
                           Species = 1,
                           Qtimes = NULL,
                           CovBinned = NULL,
                           Ncat = NULL)
{
  # TimeSim Time trajectory (vector)
  # Origin Area of origin (character)
  # CovBinned Binned covariate (vector)
  if (Origin == "random")
  {
    Start <- sample(2:3, 1)
  } else
  {
    Start <- as.numeric(Origin)
  }

  if (is.null(CovBinned))
  {
    CovBinned <- NA
  }

  SimDf <- data.frame(subject = Species, time = TimeSim, state = Start, cov = CovBinned,
                      DivA = NA_integer_, DivB = NA_integer_, Strata = 1,
                      stateSampling = NA_integer_,
                      GammaCat = ifelse(is.null(Ncat), 1, sample(1:Ncat, 1)))

  if (!is.null(Qtimes))
  {
    SimDf$Strata <- findInterval(TimeSim, Qtimes) + 1
  }
  return(SimDf)
}
