#' @title Biogeographic stochastic mapping for fossil biogeography
#'
#' @description This function performs a stochastically simulates the biogeographic history of taxa
#'
#' @param DesIn Input file as for a PyRateDES2 analysis
#' @param NumBsm Number of biogeographic stochastic maps
#' @param MaxAttempts Maximum number of attempts to stochastically simulate the biogeographic history of a single taxon
#' @param Step Time steps for the simulation
#' @param SimD Dispersal rates
#' @param SimE Extinction rates
#' @param Qtimes Shift times in dispersal, extinction, and sampling
#' @param VarD Strengths of covariate dependent dispersal
#' and in case of logistic correlation the midpoints
#' @param VarE Strengths of covariate dependent extinction
#' and in case of logistic correlation the midpoints
#' @param Covariate data.frame with the covariate for
#' e.g. environmental dependent dispersal or extinction
#' @param DivD Strengths of diversity dependent dispersal
#' @param DivE Strengths of diversity dependent extinction
#' @param DdE Strengths of dispersal dependent extinction
#' @param Cor Correlation with environment or diversity \cr Options 'exponential' or 'logistic'
#' @param DataInArea Simulate dynamic in that area
#' @param TraitD Data.frame or matrix of a trait influencing dispersal (first column taxon index and second the trait). Either a continuous trait or a binary trait codes as 1 and exp(1) (internal log transformation)
#' @param VarTraitD Strengths of trait dependent dispersal
#' @param TraitE Continuous influencing extinction (first column taxon index and second the trait). See TraitD for details
#' @param VarTraitE Strengths of trait dependent extinction
#' @param Verbose Should progress be printed?
#'
#' @return A list with stochastically mapped biogeographic histories
#'
#' @author Torsten Hauffe
#'
#' @export sim_DES
bsm_DES <- function(DesIn,
                    NumBsm = 1,
                    MaxAttempts = 1000,
                    Step,
                    SimD,
                    SimE,
                    Qtimes = NULL,
                    VarD = NULL,
                    VarE = NULL,
                    Covariate = NULL,
                    DivD = NULL,
                    DivE = NULL,
                    DdE = NULL,
                    Cor = "exponential",
                    DataInArea = NULL,
                    TraitD = NULL,
                    VarTraitD = NULL,
                    TraitE = NULL,
                    VarTraitE = NULL,
                    Verbose = FALSE)
{
  Nspecies <- nrow(DesIn)
  Time <- colnames(DesIn)[-1]
  Time <- as.numeric(gsub("X", "", Time))
  BinSize <- Time[1] - Time[2]
  SimQ <- rep(100, 2)
  BsmList <- vector(mode = "list", NumBsm)
  if (!is.null(Qtimes))
  {
    SimQ <- rep(SimQ, length(Qtimes))
  }
  for (y in 1:NumBsm)
  {
    if (!is.null(DivD) | !is.null(DivE) | !is.null(DdE))
    {
      stop("Stochastic mapping with diversity dependence is lasting too long")
    } else
    {
      Bsm <- matrix(NA_integer_, nrow = Nspecies, ncol = length(Time))
      rownames(Bsm) <- DesIn[, 1]
      colnames(Bsm) <- Time
      N <- 1
      TraitDSp <- NULL
      TraitESp <- NULL
      if (!is.null(VarTraitD) | !is.null(VarTraitE))
      {
        N <- 2
        if (!is.null(VarTraitD))
        {
          TraitD <- log(TraitD[, 2])
          TraitD <- TraitD - mean(TraitD)
        }
        if (!is.null(VarTraitE))
        {
          TraitE <- log(TraitE[, 2])
          TraitE <- TraitE - mean(TraitE)
        }
      }
      for (i in 1:Nspecies)
      {
        ObsDistr <- unlist(DesIn[i, -1])
        FirstObs <- max(which(is.na(ObsDistr))) + 1
        Observation <- matrix(c(Time[FirstObs], ObsDistr[FirstObs]), ncol = 2)
        CompatWithObs <- FALSE
        if (!is.null(VarTraitD))
        {
          TraitDSp <- cbind(1:2, c(exp(2 * TraitD[i]), 1))
          Observation <- matrix(c(rep(Time[FirstObs], 2),
                                  rep(ObsDistr[FirstObs], 2)),
                                ncol = 2, nrow = 2)
        }
        if (!is.null(VarTraitE))
        {
          TraitDSp <- cbind(1:2, c(exp(2 * TraitE[i]), 1))
          Observation <- matrix(c(rep(Time[FirstObs], 2),
                                  rep(ObsDistr[FirstObs], 2)),
                                ncol = 2, nrow = 2)
        }
        Counter <- 0
        while (!CompatWithObs & Counter < MaxAttempts)
        {
          Counter <- Counter  + 1
          SimRes <- sim_DES(Time = max(Time) - BinSize,
                            Step = Step,
                            BinSize = BinSize,
                            Nspecies = N,
                            SimD = SimD,
                            SimE = SimE,
                            SimQ = SimQ,
                            Qtimes = Qtimes,
                            Origin = 0,
                            VarD = VarD,
                            VarE = VarE,
                            Covariate = Covariate,
                            DivD = DivD,
                            DivE = DivE,
                            DdE = DdE,
                            Cor = "exponential",
                            DataInArea = DataInArea,
                            Observation = Observation,
                            GlobExt = NULL,
                            TraitD = TraitDSp,
                            VarTraitD = VarTraitD,
                            TraitE = TraitE,
                            VarTraitE = VarTraitE)
          SimDistr <- unlist(SimRes[[1]][1, -1])
          CompatWithObs <- sim_compatible_observed(SimDistr, ObsDistr)
          if (Verbose & (Counter %% 100 == 0)) {
            print(paste("Species:", i, "Attempt:", Counter))
          }
        }
        if (CompatWithObs)
        {
          Bsm[i, ] <- SimDistr
          if (Verbose)
          {
            print(paste("Species", i, "done.", Counter, "attempts needed"))
          }
        } else
        {
          if (Verbose)
          {
            print(paste("Species", i, "no success after", Counter, "attempts."))
          }
        }
      }
    }
    BsmList[[y]] <- Bsm
  }
  return(BsmList)
}
