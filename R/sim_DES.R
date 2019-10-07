#' @title Simulate fossil biogeography
#'
#' @description This function simulates dispersal, extinction, and sampling through time
#'
#' @param Time Positive number indicating the time covered by the simulation
#' @param Step Time steps for the simulation
#' @param BinSize Size of the time bins for the simulated data
#' @param Nspecies Number of species
#' @param SimD Dispersal rates
#' @param SimE Extinction rates
#' @param SimQ Sampling rates
#' @param Qtimes Shift times in dispersal, extinction, and sampling
#' @param Origin Area of origin \cr Options 'random', '1', or '2'
#' @param VarD Strength of covariate dependent dispersal
#' @param VarE Strength of covariate dependent extinction
#' @param Covariate data.frame with the covariate for
#' e.g. environmental dependent dispersal or extinction
#' @param DivD Strength of diversity dependent dispersal
#' @param DivE Strength of diversity dependent extinction
#' @param Cor Correlation with environment or diversity \cr Options 'linear' or 'exponential'
#' @param DataInArea Simulate dynamic in that area
#' @param Ncat  Number of categories for the discretized Gamma distribution
#' to simulate heterogeneous sampling
#' @param alpha i.e. the shape and rate of the Gamma distribution
#'
#' @return A data.frame
#'
#' @author Torsten Hauffe
#'
#' @export sim_DES
sim_DES <- function(Time,
                    Step,
                    BinSize,
                    Nspecies,
                    SimD,
                    SimE,
                    SimQ,
                    Qtimes = NULL,
                    Origin = "random",
                    VarD = NULL,
                    VarE = NULL,
                    Covariate = NULL,
                    DivD = NULL,
                    DivE = NULL,
                    Cor = "linear",
                    DataInArea = NULL,
                    Ncat = NULL,
                    alpha = NULL)
{
  if (Step > Time)
  {
    stop("Step size larger than time")
  }
  if (Step > BinSize)
  {
    stop("Step size larger than bin size")
  }
  if (BinSize > Time)
  {
    stop("BinSize larger than time")
  }
  if (Step < 1e-10)
  {
    stop("Step size small 1e-10 not possible")
  }
  if (!Origin %in% c("random", "1", "2"))
  {
    stop("Origin should be random, 1, or 2")
  }
  if ( (!is.null(VarD) & length(SimD) > 2) | (!is.null(DivD) & length(SimD) > 2) )
  {
    stop("Covariate/Diversity dependent dispersal and shifts in dispersal rate not compatible")
  }
  if ( (!is.null(VarE) & length(SimE) > 2) | (!is.null(DivE) & length(SimE) > 2) )
  {
    stop("Covariate/Diversity dependent extinction and shifts in extinction rate not compatible")
  }
  if ( sum(is.null(Ncat), is.null(alpha)) == 1 )
  {
    stop("Ncat and alpha need to be provided for sampling heterogeneity")
  }
  if (!is.null(DataInArea))
  {
    if (DataInArea == 1)
    {
      Origin <- "2"
    }
    if (DataInArea == 2)
    {
      Origin <- "1"
    }
  }

  if (!is.null(Qtimes))
  {
    Qtimes <- Time - Qtimes
  }
  TimeSim <- seq(0, Time, by = Step)
  Decimals <- nchar(gsub("(.*\\.)|([0]*$)", "", as.character(Step)))
  TimeSim <- round(TimeSim, digits = Decimals)
  SimDf <- gen_sim_df(TimeSim, Nspecies, Origin, Qtimes, Covariate, DataInArea, Ncat)
  SimDf <- sim_core(SimDf, SimD, SimE)
  SimDf <- sim_sampling(SimDf, SimQ, Step, Ncat, alpha)
  SimDfBinned <- bin_sim(SimDf, BinSize, TimeSim)
  DesInput <- get_DES_input(SimDfBinned, Time, BinSize, Nspecies, Distribution = "stateSampling")
  Res <- vector("list", length = 2)
  Res[[1]] <- DesInput
  SimDf$time <- Time - SimDf$time
  SimDf$state <- as.numeric(SimDf$state) - 1
  SimDf$stateSampling <- as.numeric(SimDf$stateSampling) - 1
  colnames(SimDf) <- c("Species", "Time", "RangeSim",
                       "Covariate", "DiversityA", "DiversityB",
                       "Strata", "RangeObs", "GammaCat")
  Res[[2]] <- SimDf
  return(Res)
}
