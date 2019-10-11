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
#' @param VarD Strengths of covariate dependent dispersal
#' and in case of logistic correlation the midpoints
#' @param VarE Strengths of covariate dependent extinction
#' and in case of logistic correlation the midpoints
#' @param Covariate data.frame with the covariate for
#' e.g. environmental dependent dispersal or extinction
#' @param DivD Strengths of diversity dependent dispersal
#' @param DivE Strengths of diversity dependent extinction
#' @param Cor Correlation with environment or diversity \cr Options 'linear' or 'exponential'
#' @param DataInArea Simulate dynamic in that area
#' @param Ncat  Number of categories for the discretized Gamma distribution
#' to simulate heterogeneous sampling
#' @param alpha i.e. the shape and rate of the Gamma distribution
#' @param Observation Two column matrix or data.frame of first observation time and area
#'
#' @return A list with four elements
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
                    alpha = NULL,
                    Observation = NULL)
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
  # if (!Origin %in% c("random", "1", "2"))
  # {
  #   stop("Origin should be random, 1, or 2")
  # }
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
  SimDf <- gen_sim_df(TimeSim, Nspecies, Origin, Qtimes,
                      Covariate, DataInArea, Ncat, Observation)
  SimDf <- sim_core(SimDf, SimD, SimE, VarD, VarE, DivD, DivE, Cor)
  SimDf <- sim_sampling(SimDf, SimQ, Step, Ncat, alpha, DataInArea)
  SimDfBinned <- bin_sim(SimDf, BinSize, TimeSim)
  DesInput <- get_DES_input(SimDfBinned, Time, BinSize, Nspecies,
                            Distribution = "stateSampling", DataInArea)
  Res <- vector("list", length = 2)
  Res[[1]] <- DesInput
  SimDf$time <- Time - SimDf$time
  SimDf$state <- as.numeric(SimDf$state) - 1
  SimDf$stateSampling <- as.numeric(SimDf$stateSampling) - 1
  colnames(SimDf) <- c("Species", "Time", "RangeSim",
                       "Covariate", "DiversityA", "DiversityB", "DiversityAB",
                       "Strata", "RangeObs", "GammaCat")
  Res[[2]] <- SimDf
  Rich <- unique(SimDf[, c("Time", "DiversityA", "DiversityB", "DiversityAB")])
  Rich <- Rich[order(Rich$Time, decreasing = TRUE), ]
  RichMat <- data.frame(Time = rev(TimeSim),
                        SimulatedDivA = 0,
                        SimulatedDivB = 0,
                        SimulatedDivAB = 0)
  RichMat[(nrow(RichMat) - nrow(Rich) + 1):nrow(RichMat), 2:4] <- Rich[, 2:4]
  Res[[3]] <- RichMat
  RichDesInput <- data.frame(Time = as.numeric(colnames(DesInput)[-1]),
                             ObservedDivA = apply(DesInput[,-1], 2, function(x) sum(x %in% c(1, 3))),
                             ObservedDivB = apply(DesInput[,-1], 2, function(x) sum(x %in% c(2, 3))),
                             ObservedDivAB = apply(DesInput[,-1], 2, function(x) sum(x %in% 3)),
                             row.names = NULL)
  Res[[4]] <- RichDesInput
  return(Res)
}
