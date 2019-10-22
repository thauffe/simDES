bin_covariate <- function(TimeSim,
                          Cov)
{
  # TimeSim Time trajectory (vector)
  # Cov Covariate for dispersal extinction (data.frame with time and covariate)
  # Simulation in msm goes forward in time by increase of time
  # In DES, time is given in units before present
  Cov <- Cov[order(Cov[, 1], decreasing = TRUE), ]
  MaxTimeSim <- max(TimeSim)
  if (MaxTimeSim > max(Cov[, 1]))
  {
    stop("Time of simulation exceeds the one of the covariate")
  }
  Cov <- Cov[Cov[, 1] <= MaxTimeSim, ]
  Cov[, 1] <- max(Cov[, 1]) - Cov[, 1]
  CovApp <- approx(Cov, xout = sort(unique(c(TimeSim, Cov[, 1]))), rule = 2 )
  Cov <- data.frame(Age = CovApp$x, Cov = CovApp$y)
  ToBin <- data.frame(Cov, Bin = findInterval(Cov[,1], TimeSim))
  BinnedCov <- aggregate(ToBin[, 2], by = list(ToBin$Bin), FUN = mean)$x
  BinnedCov <- BinnedCov - BinnedCov[length(BinnedCov)]
  return(BinnedCov)
}
