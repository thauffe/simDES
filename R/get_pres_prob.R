get_pres_prob <- function(SA, SB)
{
  P <- matrix(0, 4, 4)
  P[1, ] <- c(1, 0, 0, 0)
  P[2, ] <- c(SA, 1-SA, 0, 0)
  P[3, ] <- c(SB, 0, 1-SB, 0)
  P[4, ] <- c(SA * SB, SB * (1 - SA), SA * (1 - SB), (1 - SA) * (1 - SB))
  return(P)
}
