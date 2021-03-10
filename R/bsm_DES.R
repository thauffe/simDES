#' @title Biogeographic stochastic mapping for fossil biogeography
#'
#' @description This function performs a stochastically simulates the biogeographic history of taxa
#'
#' @param DesIn Input file as for a PyRateDES2 analysis
#' @param NumBsm Number of biogeographic stochastic maps
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
      for (i in 1:Nspecies)
      {
        ObsDistr <- unlist(DesIn[i, -1])
        FirstObs <- max(which(is.na(ObsDistr))) + 1
        Observation <- matrix(c(Time[FirstObs], ObsDistr[FirstObs]), ncol = 2)
        CompatWithObs <- FALSE
        Counter <- 0
        while (!CompatWithObs)
        {
          Counter <- Counter  + 1
          SimRes <- sim_DES(Time = max(Time) - BinSize,
                            Step = Step,
                            BinSize = BinSize,
                            Nspecies = 1,
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
                            TraitD = TraitD,
                            VarTraitD = VarTraitD,
                            TraitE = TraitE,
                            VarTraitE = VarTraitE)
          SimDistr <- unlist(SimRes[[1]][, -1])
          CompatWithObs <- sim_compatible_observed(SimDistr, ObsDistr)
          if (Verbose & (Counter %% 100 ==0)) {
            print(paste("Species:", i, "Attempt:", Counter))
          }
        }
        Bsm[i, ] <- SimDistr
        if (Verbose)
        {
          print(paste("Species", i, "done"))
        }
      }
    }
    BsmList[[y]] <- Bsm
  }
  return(BsmList)
}
