#' @title Number of dispersal and extinction events
#'
#' @description This function calculates for each taxon the number of dispersal and extinction events of a stochastic mapping of its biogeographic history
#'
#' @param Bsm The output of bsm_DES: A matrix (in a case of a single biogeographic stochastic map) or a list of matrices (in a case of several biogeographic stochastic maps)
#'
#' @return A list of two data.frames with the number of dispersal and extinction events per stochastic mapping
#'
#' @author Torsten Hauffe
#'
#' @export sim_DES
get_num_dis_ext_events <- function (Bsm)
{
  if (class(Bsm)[1] != "list")
  {
    Bsm <- list(Bsm)
  }
  ExtEv <- DisEv <- matrix(0, nrow = nrow(Bsm[[1]]), ncol = length(Bsm))
  rownames(ExtEv) <- rownames(DisEv) <- rownames(Bsm[[1]])
  for (i in 1:length(Bsm))
  {
    DisEv[, i] <- apply(Bsm[[i]], 1, function(x) sum(diff(x) > 0, na.rm = TRUE))
    ExtEv[, i] <- apply(Bsm[[i]], 1, function(x) sum(diff(x) < 0, na.rm = TRUE))
  }
  Events <- vector(mode = "list", length = 2)
  Events[[1]] <- DisEv
  Events[[2]] <- ExtEv
  names(Events) <- c("Dispersal_events", "Extinction_events")
  return(Events)
}

