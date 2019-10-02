bin_sim <- function(SimDf, BinSize, TimeSim)
{
  Time <- max(TimeSim)
  BinnedTime <- seq(0, Time, by = BinSize)
  BinsDesInput <- findInterval(TimeSim[-length(TimeSim)], BinnedTime)
  BinsDesInput <- c(BinsDesInput, max(BinsDesInput) + 1)
  SimDf$BinsDesInput <- BinsDesInput[match(SimDf$time, TimeSim)]
  SimDfBinned <- aggregate(SimDf[, c(3, 8)],
                           by = list(SimDf$BinsDesInput, SimDf$subject),
                           FUN = get_binned_range)
  colnames(SimDfBinned)[1:2] <- c("BinnedTime", "Species")
  return(SimDfBinned)
}


