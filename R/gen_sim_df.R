#' @title Template to simulate fossil biogeography
#'
#' @description This function generates a data.frame
#' to simulate dispersal and extinction through time.
#'
#' @param Time Positive number indicating the
#' @param Step Time steps for the simulation.
#' @param Nspecies Number of species
#' @param Origin Area of origin \cr Options '1', '2' or '3'
#' @param Covariate Covariate for e.g. environmental dependent dispersal or extinction
#'
#' @return A data.frame
#'
#' @author Torsten Hauffe
#'
#' @export gen_sim_df
gen_sim_df <- function(Time, Step, Nspecies, Origin = "random", Covariate = NULL)
{
  if (Step > Time)
  {
    stop("Step size larger than time")
  }
  TimeSim <- seq(0, Time, by = Step)
  LenTimeSim <- length(TimeSim)
  if (!is.null(Covariate))
  {
    BinnedCov <- bin_covariate(TimeSim, Covariate)
  }
  SimList <- vector("list", length = Nspecies)
  Decimals <- nchar(gsub("(.*\\.)|([0]*$)", "", as.character(Step)))
  for(i in 1:Nspecies)
  {
    AgeSpecies <- runif(1, min = 0, max = Time)
    AgeSpecies <- round(AgeSpecies, digits = Decimals)
    TimeSimSpecies <- seq(AgeSpecies, Time, by = Step)
    if (!is.null(Covariate))
    {
      L <- length(TimeSimSpecies)
      Covariate <- BinnedCov[(LenTimeSim - L):LenTimeSim]
    }
    SimList[[i]] <- gen_sim_df_cor(TimeSimSpecies, Origin, Species = i, CovBinned = Covariate)
  }
  SimDf <- do.call("rbind", SimList)
  SimDf <- SimDf[order(SimDf$time), ]
  return(SimDf)
}

