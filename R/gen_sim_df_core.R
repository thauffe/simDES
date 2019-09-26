gen_sim_df_cor <- function(TimeSim, Origin, Species = 1, CovBinned = NULL)
{
  # TimeSim Time trajectory (vector)
  # Origin Area of origin (character)
  # CovBinned Binned covariate (vector)
  if (Origin == "random")
  {
    Start <- sample(1:2, 1)
  } else
  {
    Start <- as.numeric(Origin)
  }
  if (is.null(CovBinned))
  {
    CovBinned <- NA
  }
  SimDf <- data.frame(subject = Species, time = TimeSim, start = Start, cov = CovBinned)
  return(SimDf)
}
