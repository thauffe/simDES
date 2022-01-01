bin_covariate <- function(TimeSim,
                          CovList,
                          CovNames)
{
  # TimeSim Time trajectory (vector)
  # CovList Covariate for dispersal extinction (data.frame with time and covariate)
  # Simulation in msm goes forward in time by increase of time
  # In DES, time is given in units before present
  LenCovList <- length(CovList)
  BinnedCovList <- vector(mode = "list", length = LenCovList)
  for (i in 1:LenCovList) {
    Cov <- CovList[[i]]
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
    BinnedCovList[[i]] <- BinnedCov - mean(BinnedCov)
  }
  BinnedCov <- do.call("cbind", BinnedCovList)
  BinnedCov <- BinnedCov[, rep(1:LenCovList, each = 2)]
  colnames(BinnedCov) <- paste0(rep(CovNames, LenCovList), "_", rep(1:LenCovList, each = 2))
  return(BinnedCov)
}
