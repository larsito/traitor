#' Plot decay of rank correlations.
#' 
#' Plots the decay of rank correlations between functional diversity indices of 
#' original and reduced community data, as calculated by the
#' \code{\link{reduceFD}} function.
#' @param reduced object created by \code{\link{reduceFD}}.
#' @param index one of the indices calculated by \code{\link{reduceFD}}.
#' @return A plot of decaying rank correlation with one line for each of the two
#'   reduction Scenarios. Only one plot for one functional diversity index is 
#'   produced at a time.
#' @export
plot_reduceFD <- function(reduced, index = c("FRich", "FEve", "FDis", "Rao", "CWM")){
  if (length(index) > 1) stop("Please choose one index at a time")
  rem_seq <- reduced$rem_seq
  x <- c(rev(rem_seq), 1)
  y_sc1 <- reduced$rank_cor_scenario1[, index]
  y_sc2 <- reduced$rank_cor_scenario2[, index]
  plot(x, y_sc1, type="n", main = index, xaxt = "n", 
       xlab='Proportion of abundance retained', ylab = "Rank correlation", 
       xlim =c(range(x)))
  axis(1, at=rev(x), labels=x)
  lines(x, y_sc1)
  lines(x, y_sc2, lty=2)
  legend (min(rem_seq), min(range(y_sc1)), c("Scenario 1", "Scenario 2"), lty = c(1,2), cex = 0.7, 
          y.intersp = 0.8, yjust = 0, xjust = 0)
}