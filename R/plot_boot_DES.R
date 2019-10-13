#' @title Plot bootstrapped fossil biogeography
#'
#' @description This function plots diversity trajectories
#' and their associated uncertainties through time
#'
#' @param BootDes Bootstrapped fossil biogeography from boot_DES
#' @param DivToPlot Diversities to plot
#' @param Col Colors for the different diversities
#' @param Transparency Transparency for the confidence intervals
#'
#' @return A plot
#'
#' @author Torsten Hauffe
#'
#' @export plot_boot_DES
plot_boot_DES <- function(BootDes,
                          DivToPlot = 1:3,
                          Col = c("orange", "dodgerblue", "green"),
                          Transparency = 0.5)
{
  Ci <- unlist(BootDes[[3]])
  CiNames <- names(Ci)
  Divs <- sapply(DivToPlot, function(x) max(Ci[CiNames[grepl(c("A", "B", "AB")[x], CiNames )]]))
  Ylim <- ceiling(max(Divs))
  TimeSim <- as.numeric(rownames(BootDes[[2]]))
  plot(1, 1, type = "n", ylim = c(0, Ylim), xlim = c(max(TimeSim), min(TimeSim)),
       xlab = "Time", ylab = "Diversity")
  LenCi <- length(BootDes[[3]])
  for(d in DivToPlot)
  {
    for(i in 1:LenCi)
    {
      polygon(c(TimeSim, rev(TimeSim)),
              c(BootDes[[3]][[i]][[d]][, 1], rev(BootDes[[3]][[i]][[d]][, 2])),
              col = adjustcolor(Col[d], alpha = Transparency/LenCi),
              border = NA)
    }
  }
  for(d in DivToPlot)
  {
    lines(TimeSim, BootDes[[2]][, d], col = Col[d])
  }
}
