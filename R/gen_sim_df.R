#' @title Template to simulate fossil biogeography
#'
#' @description This function generates a data.frame
#' to simulate dispersal and extinction through time.
#'
#' @param TimeSim Time steps for the simulation
#' @param Nspecies Number of species
#' @param Origin Area of origin \cr Options 'random', '1', or '2'
#' @param Covariate data.frame with the covariate for
#' e.g. environmental dependent dispersal or extinction
#'
#' @return A data.frame
#'
#' @author Torsten Hauffe
#'
#' @export gen_sim_df
gen_sim_df <- function(TimeSim,
                       Nspecies,
                       Origin = "random",
                       Qtimes = NULL,
                       Covariate = NULL,
                       DataInArea = NULL,
                       Ncat = NULL)
{
  LenTimeSim <- length(TimeSim)
  Time <- max(TimeSim)
  if (!is.null(Covariate))
  {
    BinnedCov <- bin_covariate(TimeSim, Covariate)
  }
  SimList <- vector("list", length = Nspecies)
  for(i in 1:Nspecies)
  {
    if (is.null(DataInArea))
    {
      AgeSpecies <- runif(1, min = 0, max = Time)
      Ts <- findInterval(AgeSpecies, TimeSim)
      TimeSimSpecies <- TimeSim[Ts:LenTimeSim]

    } else
    {
      TimeSimSpecies <- TimeSim
      Origin <- DataInArea + 1
    }

    if (!is.null(Covariate))
    {
      L <- length(TimeSimSpecies)
      CovBinnedSpecies <- BinnedCov[(LenTimeSim - L + 1):LenTimeSim]
    }
    else
    {
      CovBinnedSpecies <- NULL
    }

    SimList[[i]] <- gen_sim_df_cor(TimeSimSpecies, Origin, Species = i,
                                   Qtimes, CovBinned = CovBinnedSpecies, Ncat)
  }
  SimDf <- do.call("rbind", SimList)
  SimDf <- SimDf[order(SimDf$subject, SimDf$time), ]
  return(SimDf)
}

