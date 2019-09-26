#' @title Discretized Gamma distribution
#'
#' @description This function discretizes the Gamma distribution and
#' is used to simulate heterogeneous fossil preservation via multiplication
#' with the preservation rate
#'
#' @param alpha i.e. the shape and rate of the Gamma distribution
#' @param Ncat Number of categories for the discretized Gamma distribution
#'
#' @return Discretized gamma distribution
#'
#' @author Torsten Hauffe
#'
#' @references Yang, Z. (1993).
#' Maximum-likelihood estimation of phylogeny from DNA sequences
#' when substitution rates differ over sites.
#' Molecular Biology and Evolution 10(6): 1396-1401.
#'
#' @export get_gamma_rates
get_gamma_rates <- function(alpha, Ncat)
{
  SeqGamma <- seq(0, 1, length.out = Ncat + 1)
  YangGammaQuant <- (SeqGamma - SeqGamma[2]/2)[-1]
  B <- alpha
  M <- qgamma(YangGammaQuant, alpha, B)
  S <- Ncat / sum(M)
  Res <- M * S
  return(Res)
}
