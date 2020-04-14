sim_core <- function(SimDf,
                     Dis,
                     Ext,
                     VarD = NULL,
                     VarE = NULL,
                     DivD = NULL,
                     DivE = NULL,
                     DdE = NULL,
                     Cor = "exponential")
{
  # Loop is not very efficient but the only way to
  # trace richnesses for diversity-dependence
  # Round to 10th decimal because there are sometimes problems with unique()
  SimDf$time <- round(SimDf$time, 10)
  UniqueTime <- sort(unique(SimDf$time))
  IdxDiv <- SimDf$time == min(UniqueTime)
  SimDf[IdxDiv, "DivA"] <- sum(SimDf[IdxDiv, "state"] == 2)
  SimDf[IdxDiv, "DivB"] <- sum(SimDf[IdxDiv, "state"] == 3)
  SimDf[IdxDiv, "DivAB"] <- sum(SimDf[IdxDiv, "state"] == 4)
  for(i in 2:length(UniqueTime))
  {
    TimeCovered <- UniqueTime[(i - 1):i]
    Species <- SimDf[SimDf$time %in% TimeCovered, "subject"]
    # Species <- SimDf[SimDf$time == UniqueTime[i - 1], "subject"]
    # Remove species that occur just once over this time step
    # (may happen because of global extinction)
    SpeciesFreq <- table(Species)
    Species <- as.numeric(names(SpeciesFreq[SpeciesFreq > 1]))
    if (length(Species) > 0)
    {
      Start <- SimDf[SimDf$time == UniqueTime[i - 1] & SimDf$subject %in% Species, "state"]
      Strata <- SimDf[SimDf$time == UniqueTime[i], "Strata"][1]
      IdxDES <- (2 * Strata - 1):(2 * Strata)
      # No covariation with dispersal
      if (is.null(VarD) & is.null(DivD))
      {
        D <- Dis[IdxDES]
      } else # Covariation
      {
        # Covariation with environment
        if (!is.null(VarD))
        {
          CovTmp <- SimDf[SimDf$time == UniqueTime[i - 1], "cov"][1]
          if (Cor == "exponential") # Exponential covariation
          {
            D <- Dis * exp(VarD * CovTmp)
          } else # Logistic covariation
          {
            D <- Dis / (1 + exp(-VarD[1:2] * (CovTmp - VarD[3:4])))
          }
        } else
        {
          # Diversity dependent dispersal
          # (less likely colonization if there are already many taxa in the sink area)
          DivTmp <- SimDf[SimDf$time == UniqueTime[i - 1], c("DivB","DivA")][1, ]
          # DivTmp <- log1p(unlist(DivTmp))
          DivTmp <- unlist(DivTmp)
          if (Cor == "exponential") # Exponential covariation
          {
            # D <- Dis * exp(DivD * DivTmp)
            # D <- Dis + DivD * DivTmp
            D <- Dis * (1 - (DivTmp/DivD))
            D[D < 0] <- 0

          } else # Logistic covariation
          {
            D <- Dis / (1 + exp(-DivD[1:2] * (DivTmp - DivD[3:4])))
          }
        }
      }
      # No covariation with extinction
      if (is.null(VarE) & is.null(DivE) & is.null(DdE))
      {
        E <- Ext[IdxDES]
      } else # Covariation
      {
        # Covariation with environment
        if (!is.null(VarE))
        {
          CovTmp <- SimDf[SimDf$time == UniqueTime[i - 1], "cov"][1]
          if (Cor == "exponential") # Exponential covariation
          {
            E <- Ext * exp(VarE * CovTmp)
          } else # Logistic covariation
          {
            E <- Ext / (1 + exp(-VarE[1:2] * (CovTmp - VarE[3:4])))
          }
        }
        # Diversity dependent extinction
        if (!is.null(DivE))
        {
          # Diversity dependent extinction
          # (more likely extinction if there are already many taxa in the focal area)
          DivTmp <- SimDf[SimDf$time == UniqueTime[i - 1], c("DivA","DivB")][1, ]
          # DivTmp <- log1p(unlist(DivTmp))
          DivTmp <- unlist(DivTmp)
          if (Cor == "exponential") # Exponential covariation
          {
            # E <- Ext * exp(DivE * DivTmp)
            E <- Ext + DivE * DivTmp
            # E <- Ext/(1 - (DivTmp/DivE))
            # E[E < 0] <- 0
          } else # Logistic covariation
          {
            E <- Ext / (1 + exp(-DivE[1:2] * (DivTmp - DivE[3:4])))
          }
        }
        # Dispersal dependent extinction
        else
        {
          DisTmp <- SimDf[SimDf$time == UniqueTime[i - 1], c("num_d21", "num_d12")][1, ]
          DisTmp <- unlist(DisTmp)
          DivTmp <- SimDf[SimDf$time == UniqueTime[i - 1], c("DivA","DivB")][1, ]
          DivTmp <- unlist(DivTmp)
          # E <- Ext * exp(DdE * DisTmp/(DivTmp + 1))
          E <- Ext + DdE * DisTmp/(DivTmp + 1)
        }

      }
      Q <- make_Q(D, E)
      Index <- SimDf$time %in% TimeCovered & SimDf$subject %in% Species
      StatesBeforeTransition <- SimDf[Index, "state"]
      SimDf[Index, "state"] <- simmulti.msm(SimDf[Index, ],
                                            qmatrix = Q,
                                            start = Start)$state
      IdxDiv <- SimDf$time == UniqueTime[i]
      SimDf[IdxDiv, "DivA"] <- sum(SimDf[IdxDiv, "state"] %in% c(2, 4))
      SimDf[IdxDiv, "DivB"] <- sum(SimDf[IdxDiv, "state"] %in% c(3, 4))
      SimDf[IdxDiv, "DivAB"] <- sum(SimDf[IdxDiv, "state"] %in% 4)
      IdxRate <- SimDf$time == UniqueTime[i - 1]
      RepIdxRate <- sum(IdxRate)
      SimDf[IdxRate, c("rate_d12", "rate_d21")] <- rep(D, each = RepIdxRate)
      SimDf[IdxRate, c("rate_e1", "rate_e2")] <- rep(E, each = RepIdxRate)
      IdxDis <- SimDf$time == UniqueTime[i]
      SimDf[IdxDis, "num_d12"] <- sum(SimDf[Index, "state"] - StatesBeforeTransition == 2)
      SimDf[IdxDis, "num_d21"] <- sum(SimDf[Index, "state"] - StatesBeforeTransition == 1)
    }
  }
  return(SimDf)
}

