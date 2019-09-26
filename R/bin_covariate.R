bin_covariate <- function(TimeSim, Cov)
{
  # TimeSim Time trajectory (vector)
  # Cov Covariate for dispersal extinction (data.frame with time and covariate)
  # Simulation in msm goes forward in time by increase of time
  # In DES, time is given in units before present
  Cov <- Cov[nrow(Cov):1, ]
  Cov[, 1] <- Cov[, 1] - max(Cov[, 1])
  MaxTimeSim <- max(TimeSim)
  if (MaxTimeSim > max(Cov[, 1]))
  {
    stop("Time of simulation exceeds the one of the covariate")
  }
  Cov <- Cov[Cov[,1] <= MaxTimeSim, ]
  ToBin <- data.frame(Cov, Bin = findInterval(Cov[,1], TimeSim))
  BinnedCov <- aggregate(ToBin[, 2], by = list(ToBin$Bin), FUN = mean)$x
  return(BinnedCov)
}
