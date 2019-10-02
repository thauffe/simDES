make_Q <- function(D, E)
{
  Q <- rbind(c(0,     0,     0,     0),
             c(E[1], 0,     0,     D[1]),
             c(E[2], 0,     0,     D[2]),
             c(0,     E[2], E[1], 0))
  return(Q)
}
