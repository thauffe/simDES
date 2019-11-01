#' @title Bootstrap fossil biogeography
#'
#' @description This function simulates several times the biogeographic events
#' dispersal, extinction, and sampling to obtain diversity trajectories
#' and their associated uncertainties through time
#'
#' @param Nsim Number of bootstrap simulations
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
#' @param Cor Correlation with environment or diversity \cr Options 'exponential' or 'logistic'
#' @param DataInArea Simulate dynamic in that area
#' @param Ncat  Number of categories for the discretized Gamma distribution
#' to simulate heterogeneous sampling
#' @param alpha i.e. the shape and rate of the Gamma distribution
#' @param Observation Two column matrix or data.frame of first observation time and area
#' @param ConfInt Confidence interval for simualted diversities
#' (e.g. c(0.95, 0.68) for two and one standard deviation)
#' @param Ncores Number of CPU cores to use for simulations
#'
#' @return A list of three elements.
#'
#' @author Torsten Hauffe
#'
#' @export boot_DES
boot_DES <- function(Nsim = 10,
                     Time,
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
                     Cor = "exponential",
                     DataInArea = NULL,
                     Ncat = NULL,
                     alpha = NULL,
                     Observation = NULL,
                     ConfInt = 0.95,
                     Ncores = 1)
{
  TimeSim <- rev(seq(0, Time, by = Step))
  DivA <- DivB <- DivAB <- matrix(NA_integer_,
                                  nrow = length(TimeSim),
                                  ncol = Nsim)
  rownames(DivA) <- rownames(DivB) <- rownames(DivAB) <- TimeSim
  Res <- vector("list", length = 3)
  registerDoParallel(Ncores)
  SimTmp <- foreach(iter = 1:Nsim,
                    .inorder = FALSE) %dopar% sim_DES(Time,
                                                      Step,
                                                      BinSize,
                                                      Nspecies,
                                                      SimD,
                                                      SimE,
                                                      SimQ,
                                                      Qtimes,
                                                      Origin,
                                                      VarD,
                                                      VarE,
                                                      Covariate,
                                                      DivD,
                                                      DivE,
                                                      Cor,
                                                      DataInArea,
                                                      Ncat,
                                                      alpha,
                                                      Observation)[[3]]
  stopImplicitCluster()
  for(i in 1:Nsim)
  {
      DivA[, i] <- SimTmp[[i]]$SimulatedDivA
      DivB[, i] <- SimTmp[[i]]$SimulatedDivB
      DivAB[, i] <- SimTmp[[i]]$SimulatedDivAB
  }
  DivList <- vector("list", length = 3)
  DivList[[1]] <- DivA
  DivList[[2]] <- DivB
  DivList[[3]] <- DivAB
  names(DivList) <- c("SimulatedDivA", "SimulatedDivB", "SimulatedDivAB")
  Res[[1]] <- DivList

  DivMean <- data.frame(MeanSimDivA = rowMeans(DivA),
                        MeanSimDivB = rowMeans(DivB),
                        MeanSimDivAB = rowMeans(DivAB))
  rownames(DivMean) <- TimeSim
  Res[[2]] <- DivMean

  ConfInt <- sort(ConfInt, decreasing = TRUE)
  Probs1 <- (1 - ConfInt)/2
  Probs2 <- Probs1 + ConfInt
  Probs <- cbind(Probs1, Probs2)
  CiList <- vector("list", length = length(ConfInt))
  for(i in 1:length(ConfInt)){
    Ci <- vector("list", length = 3)
    Ci[[1]] <- t(apply(DivA, 1, function(x) quantile(x, Probs[i, ])))
    Ci[[2]] <- t(apply(DivB, 1, function(x) quantile(x, Probs[i, ])))
    Ci[[3]] <- t(apply(DivAB, 1, function(x) quantile(x, Probs[i, ])))
    rownames(Ci[[1]]) <- rownames(Ci[[2]]) <- rownames(Ci[[3]]) <- TimeSim
    names(Ci) <- paste(c("A", "B", "AB"), ConfInt[i], sep = "_")
    CiList[[i]] <- Ci
  }
  Res[[3]] <- CiList
  return(Res)
}
