sim_core <- function(SimDf,
                     Dis,
                     Ext,
                     VarD = NULL,
                     VarE = NULL,
                     DivD = NULL,
                     DivE = NULL,
                     Cor = "linear")
{
  # Loop is not very efficient but the only way to
  # trace richnesses for diversity-dependence
  # Round to 10th decimal because there are sometimes problems with unique()
  SimDf$time <- round(SimDf$time, 10)
  UniqueTime <- sort(unique(SimDf$time))
  IdxDiv <- SimDf$time == min(UniqueTime)
  SimDf[IdxDiv, "DivA"] <- sum(SimDf[IdxDiv, "state"] == 2)
  SimDf[IdxDiv, "DivB"] <- sum(SimDf[IdxDiv, "state"] == 3)
  for(i in 2:length(UniqueTime))
  {
    TimeCovered <- UniqueTime[(i - 1):i]
    Species <- SimDf[SimDf$time == UniqueTime[i - 1], "subject"]
    Start <- SimDf[SimDf$time == UniqueTime[i - 1] & SimDf$subject %in% Species, "state"]
    Strata <- SimDf[SimDf$time == UniqueTime[i], "Strata"][1]
    IdxDES <- (2 * Strata - 1):(2 * Strata)
    if (is.null(VarD))
    {
      D <- Dis[IdxDES]
    } else
    {
      D <- Dis[1:2]
    }
    if (is.null(VarE))
    {
      E <- Ext[IdxDES]
    } else
    {
      E <- Ext[1:2]
    }
    Q <- make_Q(D, E)
    Index <- SimDf$time %in% TimeCovered & SimDf$subject %in% Species
    SimDf[Index, "state"] <- simmulti.msm(SimDf[Index, ],
                                          qmatrix = Q,
                                          start = Start)$state
    IdxDiv <- SimDf$time == UniqueTime[i]
    SimDf[IdxDiv, "DivA"] <- sum(SimDf[IdxDiv, "state"] %in% c(2, 4))
    SimDf[IdxDiv, "DivB"] <- sum(SimDf[IdxDiv, "state"] %in% c(3, 4))
  }
  return(SimDf)
}
