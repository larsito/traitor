#' Functional diversity indices for original and reduced community data.
#' 
#' \code{reduceFD} uses \code{\link[FD]{dbFD}} to calculate functional diversity
#' indices for original and reduced community data for two different scenarios 
#' (pool-wise and plot-wise). Also, rank correlations between the FD indices of 
#' original and reduced data are calculated (Majekova et al., unpubl.).
#' @param com a matrix containing the community data, with species as columns 
#'   and plots (samples) as rows.
#' @param traits a vector, data frame or matrix containing trait values for the 
#'   species in \code{com}.
#' @param reduction the relative abundance removed in each reduction step.
#' @param remain the amount of total relative abundance retained in the 
#'   community. Gives an error if this value would lead to plots having less 
#'   than two species.
#' @return A list with five elements. See details.
#' @details The function returns a list of five elements, containing the 
#'   following objects: \describe{\item{\code{$rem_seq}}{A vector containing the 
#'   sequence of reduction steps used to calculate functional diversities. The 
#'   last (smallest) value in this vector is equal to \code{remain}.} 
#'   \item{\code{$FD_reduced_scenario1}}{This object contains a list of the same 
#'   length as \code{$rem_seq}, where each element contains the functional 
#'   diversity indices for each of the reduction 
#'   steps.}\item{\code{$FD_reduced_scenario1}}{The same as the previous object, 
#'   but for Scenario 2}\item{\code{$rank_cor_scenario1}}{The rank correlations 
#'   between the functional diversity indices of each of the reduction steps 
#'   (for Scenario 1) and the functional diversity indices of the original data.
#'   This value demonstrates, for each functional diversity index separately, 
#'   how the index from increasingly reduced communities contain less and less 
#'   "true" information about the actual functional diversities across the plots
#'   (samples).}\item{\code{$rank_cor_scenario2}}{The same as the previous 
#'   object, but for Scenario 2.}}
#' @seealso \code{\link[FD]{dbFD}} for details on functional diversity 
#'   calculation, \code{\link{plot.reduceFD}} for a graphical visualisation of 
#'   the output of reduceFD.
#' @references Majekova, Maria; et al. Evaluating functional diversity: missing 
#'   trait data and the importance of species abundance structure and data 
#'   transformation. Unpublished.
#' @examples
#' data(ohrazeni)
#' w <- reduceFD(ohrazeni$community[, -49], 
#'      na.omit(ohrazeni$traits[, "ch", drop = F]), 
#'      reduction = 0.1, remain = 0.7)   
#' plot_reduceFD(w)
#' @export
reduceFD <- function(com, traits, reduction = 0.05, remain = 0.5) {
  if(any(colnames(com) != rownames(traits))) {
    stop(cat("Species in community and trait data have to be the same and in the same order"))
  }  
  if(max(colSums(com)/sum(com)) > remain) {
    stop(cat("Please chose a smaller value for the remain argument"))
  }
  rem_seq <- seq(1-reduction, remain, -reduction)
  red_L <- mapply(reducer, list(com), rem_seq, SIMPLIFY = FALSE)
  FDred <- lapply(red_L, function(a) lapply(a, functDiv, traits))
  FD_org <- lapply(lapply(FDred, "[", 1), unlist, recursive = FALSE)[[1]]
  FD_org <- FD_org[c(1, 3, 5:8)]
  FD_red_sc1 <- lapply(lapply(FDred, "[", 2), unlist, recursive = FALSE)
  FD_red_sc2 <- lapply(lapply(FDred, "[", 3), unlist, recursive = FALSE)
  names(FD_red_sc1) <- names(FD_red_sc2) <- paste("remain", rem_seq, 
                                                  sep = "_")
  FD_red_sc1 <- lapply(lapply(FD_red_sc1, function(x) do.call(cbind, x)), 
                       function(x) x[, c(1,3,5:8)])
  FD_red_sc2 <- lapply(lapply(FD_red_sc2, function(x) do.call(cbind, x)), 
                       function(x) x[, c(1,3,5:8)])
  FD_red_sc1_2 <- list()
  for(i in 1:6){
    FD_red_sc1_2[[i]] <- cbind(FD_org[[i]], do.call(cbind, lapply(FD_red_sc1, "[", i)))
  }
  FD_red_sc2_2 <- list()
  for(i in 1:6){
    FD_red_sc2_2[[i]] <- cbind(FD_org[[i]], do.call(cbind, lapply(FD_red_sc2, "[", i)))
  }
  
  rank_cor_sc1 <- do.call(cbind, lapply(lapply(FD_red_sc1_2, cor, 
                                               method = "spearman", 
                                               use = "pairwise.complete.obs"), 
                                        "[", , 1))
  dimnames(rank_cor_sc1) <- list(c(1, rem_seq), c("NrSpecies", "FRich", "FEve", 
                                                  "FDis", "Rao", "CWM"))
  rank_cor_sc2 <- do.call(cbind, lapply(lapply(FD_red_sc2_2, cor, 
                                               method = "spearman", 
                                               use = "pairwise.complete.obs"), 
                                        "[", , 1))
  dimnames(rank_cor_sc2) <- list(c(1, rem_seq), c("NrSpecies", "FRich", "FEve", 
                                                  "FDis", "Rao", "CWM"))
  out <- list(rem_seq = rem_seq, 
              FD_reduced_scenario1 = FD_red_sc1, 
              FD_reduced_scenario2 = FD_red_sc2, 
              rank_cor_scenario1 = rank_cor_sc1, 
              rank_cor_scenario2 = rank_cor_sc2)
  return(out)
}