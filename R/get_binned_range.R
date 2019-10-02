get_binned_range <- function(X)
{
  if (4 %in% X | (2 %in% X & 3 %in% X))
  {
    A <- 4
  } else if (2 %in% X & ! 3 %in% X)
  {
    A <- 2
  } else if (3 %in% X)
  {
    A <- 3
  } else
  {
    A <- 1
  }
  return(A)
}
