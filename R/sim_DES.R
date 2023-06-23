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
#' @param Origin Vector with the possible areas of origin \cr Options: 0 (i.e. random, default) or any combination of 1, 2, and/or 3
#' @param VarD Strengths of covariate dependent dispersal
#' and in case of logistic correlation the midpoints
#' @param VarE Strengths of covariate dependent extinction
#' and in case of logistic correlation the midpoints
#' @param CovariateDisp List of a single or several data.frames with the covariate for
#' e.g. environmental dependent dispersal
#' @param CovariateExt List of a single or several data.frames with the covariate for
#' e.g. environmental dependent extinction
#' @param DivD Strengths of diversity dependent dispersal
#' @param DivE Strengths of diversity dependent extinction
#' @param DdE Strengths of dispersal dependent extinction
#' @param Cor Correlation with environment or diversity \cr Options 'exponential' or 'logistic'
#' @param DataInArea Simulate dynamic in that area
#' @param Ncat  Number of categories for the discretized Gamma distribution
#' to simulate heterogeneous sampling or Inf for a random draw of the Gamma distribution
#' @param alpha i.e. the shape and rate of the Gamma distribution
#' @param Observation Two column matrix or data.frame of first observation time and area
#' @param GlobExt Global extinction rate
#' @param TraitD Data.frame or matrix of a trait influencing dispersal (first column taxon index and second the trait). Either a continuous trait or a binary trait codes as 1 and exp(1) (internal log transformation)
#' @param VarTraitD Strengths of trait dependent dispersal
#' @param TraitE Continuous influencing extinction (first column taxon index and second the trait). See TraitD for details
#' @param VarTraitE Strengths of trait dependent extinction
#' @param TraitS Continuous influencing sampling (first column taxon index and second the trait). See TraitD for details
#' @param VarTraitS Strengths of trait dependent sampling
#' @param CatTraitD Data.frame or matrix of the effect of a discrete trait on dispersal (first column taxon index and second the state-dependent multiplier of baseline rate; requires the same order as Observation).
#' @param CatTraitE Data.frame or matrix of the effect of a discrete trait on extinction (first column taxon index and second the state-dependent multiplier of baseline rate; requires the same order as Observation).
#' @param CatTraitS Data.frame or matrix of the effect of a discrete trait on sampling (first column taxon index and second the state-dependent multiplier of baseline rate; requires the same order as Observation).
#' @param Verbose Should messages be printed?
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
                    Origin = 0,
                    VarD = NULL,
                    VarE = NULL,
                    CovariateDis = NULL,
                    CovariateExt = NULL,
                    DivD = NULL,
                    DivE = NULL,
                    DdE = NULL,
                    Cor = "exponential",
                    DataInArea = NULL,
                    Ncat = NULL,
                    alpha = NULL,
                    Observation = NULL,
                    GlobExt = NULL,
                    TraitD = NULL,
                    VarTraitD = NULL,
                    TraitE = NULL,
                    VarTraitE = NULL,
                    TraitS = NULL,
                    VarTraitS = NULL,
                    CatTraitD = NULL,
                    CatTraitE = NULL,
                    CatTraitS = NULL,
                    Verbose = FALSE)
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
  TimeTmp <- seq(0, Time, by = BinSize)
  if (max(TimeTmp) != Time)
  {
    BinSize = Time/length(TimeTmp)
    if (Verbose)
    {
      cat(paste("Time not the product of time bin size.\nReset to", BinSize))
    }

  }
  if (Step < 1e-10)
  {
    stop("Step size small 1e-10 not possible")
  }
  if ( !all(Origin %in% 0:3) | (any(Origin %in% 0) & any(Origin %in% 1:3)) )
  {
    stop("Origin should be 0 (i.e. random) or any combination of 1, 2, and/or 3")
  }
  if (Origin == 0)
  {
    Origin <- 1:3
  }
  if ( (!is.null(VarD) & length(SimD) > 2) | (!is.null(DivD) & length(SimD) > 2) )
  {
    stop("Covariate/Diversity dependent dispersal and shifts in dispersal rate not compatible")
  }
  if ( (!is.null(VarE) & length(SimE) > 2) | (!is.null(DivE) & length(SimE) > 2) | (!is.null(DdE) & length(SimE) > 2) )
  {
    stop("Covariate/Diversity/Dispersal dependent extinction and shifts in extinction rate not compatible")
  }
  if ( sum(is.null(Ncat), is.null(alpha)) == 1 )
  {
    stop("Ncat and alpha need to be provided for sampling heterogeneity")
  }
  if (!is.null(DataInArea) & is.null(Observation))
  {
    if (DataInArea == 1)
    {
      Origin <- 2
    }
    if (DataInArea == 2)
    {
      Origin <- 1
    }
  }

  if (!is.null(Qtimes))
  {
    Qtimes <- sort(Qtimes, decreasing = TRUE)
    Qtimes <- Time - Qtimes
  }
  TimeSim <- seq(0, Time, by = Step)
  Decimals <- nchar(gsub("(.*\\.)|([0]*$)", "", as.character(Step)))
  TimeSim <- round(TimeSim, digits = Decimals)
  SimDf <- gen_sim_df(TimeSim, Nspecies, Origin, Qtimes,
                      CovariateDis, CovariateExt, DataInArea, Ncat, Observation, GlobExt)
  SimDf <- sim_core2(SimDf, SimD, SimE,
                     VarD, VarE, DivD, DivE, DdE, Cor,
                     TraitD, VarTraitD, TraitE, VarTraitE, CatTraitD, CatTraitE)
  SimDf <- sim_sampling(SimDf, SimQ, Nspecies, Step, Ncat, alpha, DataInArea, TraitS, VarTraitS, CatTraitS)
  SimDfBinned <- bin_sim(SimDf, BinSize, TimeSim)
  DesInput <- get_DES_input(SimDfBinned, Time, BinSize,
                            Distribution = "stateSampling", DataInArea)
  Res <- vector("list", length = 2)
  Res[[1]] <- DesInput
  SimDf$time <- Time - SimDf$time
  SimDf$state <- as.numeric(SimDf$state) - 1
  SimDf$stateSampling <- as.numeric(SimDf$stateSampling) - 1
  colnames(SimDf)[1:3] <- c("Species", "Time", "RangeSim")
  L <- ncol(SimDf)
  colnames(SimDf)[(L - 12):L] <- c("DiversityA", "DiversityB", "DiversityAB",
                                   "Strata", "RangeObs", "GammaCat",
                                   "rate_d12", "rate_d21", "rate_e1", "rate_e2",
                                   "num_d12", "num_d21", "StateObserved")
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
  if ( !is.null(VarTraitD) )
  {
    TraitD[, 1] <- paste0("SP_", TraitD[, 1])
    TraitD <- TraitD[TraitD[, 1] %in% DesInput$scientificName, ]
    Res[[length(Res) + 1]] <- TraitD
  }
  if ( !is.null(VarTraitE) )
  {
    TraitE[, 1] <- paste0("SP_", TraitE[, 1])
    TraitE <- TraitE[TraitE[, 1] %in% DesInput$scientificName, ]
    Res[[length(Res) + 1]] <- TraitE
  }
  if ( !is.null(VarTraitS) )
  {
    TraitS[, 1] <- paste0("SP_", TraitS[, 1])
    TraitS <- TraitS[TraitS[, 1] %in% DesInput$scientificName, ]
    Res[[length(Res) + 1]] <- TraitS
  }
  if ( !is.null(CatTraitD) )
  {
    CatTraitD[, 1] <- paste0("SP_", CatTraitD[, 1])
    CatTraitD <- CatTraitD[CatTraitD[, 1] %in% DesInput$scientificName, ]
    Res[[length(Res) + 1]] <- CatTraitD
  }
  if ( !is.null(CatTraitE) )
  {
    CatTraitE[, 1] <- paste0("SP_", CatTraitE[, 1])
    CatTraitE <- CatTraitE[CatTraitE[, 1] %in% DesInput$scientificName, ]
    Res[[length(Res) + 1]] <- CatTraitE
  }
  return(Res)
}
