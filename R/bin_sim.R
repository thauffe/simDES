bin_sim <- function(SimDf, BinSize, TimeSim)
{
  Time <- max(TimeSim)
  BinnedTime <- seq(0, Time, by = BinSize)
  Decimals <- nchar(gsub("(.*\\.)|([0]*$)", "", as.character(BinSize)))
  BinnedTime <- round(BinnedTime, digits = Decimals)
  BinsDesInput <- findInterval(TimeSim[-length(TimeSim)], BinnedTime)
  BinsDesInput <- c(BinsDesInput, max(BinsDesInput) + 1)
  SimDf$BinsDesInput <- BinsDesInput[match(SimDf$time, TimeSim)]
  SimDfBinned <- aggregate(SimDf[, c("state", "stateSampling")],
                           by = list(SimDf$BinsDesInput, SimDf$subject),
                           FUN = get_binned_range)
  colnames(SimDfBinned)[1:2] <- c("BinnedTimeIndex", "Species")
  return(SimDfBinned)
}


# bin_sim <- function(SimDf, BinSize, TimeSim)
# {
#   Time <- max(TimeSim)
#   BinnedTime <- seq(0, Time, by = BinSize)
#   nTaxa <- max(SimDf$subject)
#   SimDfBinned <- c()
#   for (j in 1:nTaxa)
#   {
#     samples_taxon_j <- SimDf[SimDf$subject == j, c(2, 3, 9)]
#     for (i in 2:(length(BinnedTime)))
#     {
#       t1 <- BinnedTime[i - 1]
#       t2 <- BinnedTime[i]
#       ind <- samples_taxon_j$time >= t1 & samples_taxon_j$time < t2
#       samples_taxon_j_time_t <- samples_taxon_j[ind, ]
#       if (nrow(samples_taxon_j_time_t) > 0)
#       {
#         Append <- c(i - 1, j,
#                     get_binned_range(samples_taxon_j_time_t$state),
#                     get_binned_range(samples_taxon_j_time_t$stateSampling))
#         SimDfBinned <- rbind(SimDfBinned, Append)
#       }
#     }
#     Append <- c(length(BinnedTime), j,
#                 samples_taxon_j$state[length(samples_taxon_j$state)],
#                 samples_taxon_j$stateSampling[length(samples_taxon_j$stateSampling)])
#     SimDfBinned <- rbind(SimDfBinned, Append)
#   }
#   colnames(SimDfBinned)[1:2] <- c("BinnedTimeIndex", "Species",
#                                   "state", "stateSampling")
#   return(SimDfBinned)
# }
