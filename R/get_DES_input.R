get_DES_input <- function(SimDfBinned, Time, BinSize, Nspecies, Distribution = "state")
{
  BinnedTime <- seq(0, Time, by = BinSize)
  DesInput <-  data.frame(scientificName = paste0("SP_", 1:Nspecies),
                          as.data.frame(matrix(NA_integer_,
                                               nrow = Nspecies,
                                               ncol = length(BinnedTime))))
  colnames(DesInput)[-1] <- rev(BinnedTime)
  for(i in 1:Nspecies){
    Species <- SimDfBinned[SimDfBinned$Species == i, ]
    #DesInput[i, 1 + findInterval(Species$BinnedTime, BinnedTime)] <- Species[, Distribution] - 1 # 1:4 -> 0:3
    DesInput[i, 1 + Species$BinnedTime] <- Species[, Distribution] - 1 # 1:4 -> 0:3
  }
  DesInput <- DesInput[rowSums(DesInput[, -1], na.rm = TRUE) > 0, ]
  return(DesInput)
}
