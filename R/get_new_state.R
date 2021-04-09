get_new_state <- function(State, Q, dT)
{
  # print(paste("Old state:", State, "Rates:", paste(Q[State, ], collapse = " ")))
  if (State != 1)
  {
    QRate <- Q[State, ]
    RanTime <- rexp(1, rate = sum(QRate))
    if (RanTime <= dT)
    {
      State <- sample(1:4, 1, prob = QRate)
    }
  }
  # print(paste("New state:", State))
  return(State)
}
