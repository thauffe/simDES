sim_compatible_observed <- function (SimDistr, ObsDistr)
{
  Compatible <- TRUE
  Incomp0 <- any(SimDistr %in% 0 & ObsDistr %in% 1:3)
  Incomp1 <- any(SimDistr %in% 1 & ObsDistr %in% c(0, 2, 3))
  Incomp2 <- any(SimDistr %in% 2 & ObsDistr %in% c(0, 1, 3))
  LenDistr <- length(SimDistr)
  IncompPresent <- SimDistr[LenDistr] == ObsDistr[LenDistr]
  if (Incomp0 | Incomp1 | Incomp2 | IncompPresent)
  {
    Compatible <- FALSE
  }
  return(Compatible)
}
