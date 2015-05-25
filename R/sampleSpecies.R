#' Evaluate trait data availability for community data.
#' 
#' This function evaluates for what proportion of species trait data is 
#' available in a set of community data. Based on this information it also 
#' provides the user with a list of species that should be sampled in order to 
#' reach a user-defined threshold of data completeness, based on the abundance 
#' data (e.g. percentage cover, individual counts, biomass) of the species.
#' @param com a matrix containing the community data, with species as columns 
#'   and plots (samples) as rows.
#' @param avail_sp a binary vector of 0 and 1, either named with the same 
#'   species as in \code{com} or with the same number and order of species as in
#'   \code{com}. 0 indicates that no trait data is available for a species. 
#'   Alternatively, a categorical or numeric vector of trait values but 
#'   containing NAs, where NA indicates species with missing trait data.
#' @param thresh a threshold value defining the minimum amount of relative 
#'   abundance that trait data should be available for.
#' @param sequential if \code{TRUE}, species to be sampled are first determined 
#'   by Scenario 1 and any species to be sampled in addition are determined by 
#'   Scenario 2. When \code{FALSE}, Scenario 2 gives the species to be sampled 
#'   in each plot to reach the threshold in each plot. See also details.
#' @references Majekova, M., et al. Evaluating functional diversity: missing 
#'   trait data and the importance of species abundance structure and data 
#'   transformation. Unpublished.
#' @return a list with two or more elements, depending on the input data and 
#'   arguments. See details.
#' @details There are different ways in which the data set is evaluated, 
#'   according to the given data and the value of the \code{sequential} 
#'   argument. In general, two different scenarios are applied to the data. In 
#'   Scenario 1 (pool-wise scenario), the calculations are carried out on the 
#'   basis of the pooled abundances of the species across the plots. In Scenario
#'   2 (plot-wise scenario), the calculations are carried out on each of the 
#'   plots in the community data. Assessing the proportion of species that has
#'   trait data available with Scenario 1 therefor returns on single value,
#'   representing the total proportion of species' abundance for which trait
#'   data is available. With Scenario 2 this value is calculated for each of the
#'   plots.
#'   
#'   When \code{sequential} is set to \code{TRUE} the function first determines 
#'   which species need to be sampled to reach the threshold according to 
#'   Scenario 1, i.e. which species need to be sampled to ensure that the 
#'   \strong{overall} (pooled) abundance covered is at least the value set by 
#'   \code{thresh}. However, since this leads to the individual plots not 
#'   necessarily reaching this threshold, subsequently species are determined by
#'   Scenario 2 (i.e. plot-wise) that need to be sampled \strong{in addition} 
#'   to the ones determined by Scenario 1, in order to reach the threshold for 
#'   each of the plots. Such a sampling strategy would be appropriate if one 
#'   needs community trait values, i.e. it does not matter which individual plot
#'   the trait measurements come from. It can also occur that the threshold is 
#'   already reached for each plot with Scenario 1, in which case a warning is 
#'   shown and the value for Scenario 2 in the output is NULL.
#'   
#'   In the second case, \code{sequential} would be set to \code{FALSE} and then
#'   species to be sampled would be determined independent for Scenario 1 and 
#'   Scenario 2. For Scenario 1, species determined would be the same, but for 
#'   Scenario 2 species would be determined without taking into account species 
#'   already suggested by Scenario 1. Hence, the function provides a full set of 
#'   species that need to be sampled in each individual plot of the community. 
#'   This would be an appropriate sampling strategy when for instance one would 
#'   want to have plot-wise trait measurements, i.e. one could than take into 
#'   account trait variation on the between-plot level, or even intraspecific 
#'   level when several individual per species are measured per plot.
#'   
#'   All these examples discussed so far assume that species in the community 
#'   data are represented by some sort of (semi)quantitative abundance measure 
#'   (e.g. cover values, individual counts, biomass). Note that Scenario 2 can 
#'   not be calculated when the community data is presence-absence only.
#'   
#'   \code{sampleSpecies} returns a list with three elements: \describe{ 
#'   \item{\code{$Original}}{Parameters of the original community data} 
#'   \item{\code{$Scenario1}}{Parameters of the community data after applying 
#'   Scenario 1} \item{\code{$Scenario2}}{Parameters of the community data after
#'   applying Scenario 2} } Each of these objects is itself a list, containing 
#'   four objects, except for \code{$Original}, which only contains the latter 
#'   three objects: \describe{\item{\code{$sample_species}}{The species to be 
#'   sampled to reach the threshold. When \code{sequential} is \code{TRUE}, this
#'   gives a single vector of species to be sampled for both Scenario 1 and 2, 
#'   with the vector under Scenario 2 given the species to be sampled in 
#'   addition to these under Scenario 1. When \code{sequential} is \code{FALSE},
#'   this gives a list of vectors of species to be sampled, where each vector 
#'   contains the species to be sampled in a given plot. \code{NULL} implies 
#'   that the abundance threshold is being reached without sampling any 
#'   (further) species than the ones that already have trait data available.} 
#'   \item{\code{$com_available}}{The species available in the community data. 
#'   If "NA", this species is not available.} \item{\code{$available_pool}}{The 
#'   total proportional abundance available on the overall (pooled) community 
#'   level.} \item{\code{$available_plot}}{The total proportional abundance 
#'   available in each plot.} }
#' @examples
#' # example community data and vector indicating trait availability for species
#' com <- matrix(c(50, 40, 30, 20, 0, 30, 0, 10, 0, 5, 10, 0, 0, 0, 5), 3, 5)
#' colnames(com) <- as.character(c(1:5))
#' avail <- c(1, 0, 1, 0, 1)
#' 
#' # run example with different argument settings
#' sampleSpecies(com, avail, 0.95, sequential = TRUE)
#' sampleSpecies(com, avail, 0.8, sequential = FALSE)
#' 
#' # run with the ohrazeni wet meadow data set
#' data(ohrazeni)
#' sampleSpecies(ohrazeni$vegetationdata, rbinom(61, 1, 0.5), 0.95, sequential = FALSE)
#' 
#' @export

sampleSpecies <- function(com, avail_sp, thresh = 0.8, sequential = TRUE) {
  out <- vector("list", 3)
  names(out) <- c("Original", "Scenario1", "Scenario2")
  if(is.data.frame(avail_sp)) {
    avail_sp <- setNames(unlist(avail_sp), rownames(avail_sp))
  }
  if(any(is.na(avail_sp))) {
    avail_sp[which(!is.na(avail_sp))] <- 1
    avail_sp[which(is.na(avail_sp))] <- 0
  }
  if (is.null(names(avail_sp))) {
    names(avail_sp) <- colnames(com)
    warning("avail_sp is without names. Species are assumed to be in the same order as in community data.", call. = FALSE)
  }
  com_av <- com[ , as.logical(avail_sp), drop = FALSE]
  relab_com <- com/rowSums(com)
  relab_com_av <- relab_com[ , as.logical(avail_sp), drop = FALSE]
  pool <- colSums(com)
  relab_pool <- pool/sum(pool)
  relab_pool_av <- relab_pool[as.logical(avail_sp)]
  tot_pool_av <- sum(relab_pool_av)
  tot_com_av <- rowSums(as.matrix(relab_com_av))
  out$Original$com_available <- com_av
  out$Original$available_pool <- tot_pool_av 
  out$Original$available_plot <- tot_com_av
  
  # Scenario 1
  to_sample1 <- sampler(relab_pool, avail_sp, thresh)
  avail_sp1 <- c(colnames(com_av), to_sample1)
  unavail_sp1 <- setdiff(colnames(com), avail_sp1)
  com_av1 <- com
  com_av1[, unavail_sp1] <- NA
  relab_com_av1 <- relab_com[, avail_sp1, drop = FALSE]
  tot_pool_av1 <- sum(relab_pool[avail_sp1])
  tot_com_av1 <- rowSums(relab_com_av1)
  if (is.null(to_sample1)) to_sample1 <- "none"
  out$Scenario1$sample_species <- sort(to_sample1)
  out$Scenario1$com_available <- com_av1
  out$Scenario1$available_pool <- tot_pool_av1
  out$Scenario1$available_plot <- tot_com_av1
  
  # Scenario 2
  matrix_is_01 <- sum(com) == sum(replace(com, com > 0, 1))
  if (matrix_is_01) {
    warning("\nScenario 2 cannot be computed as there is no ranking within a plot, species addition is based only on the frequency across all plots\n\n")
    out[[3]] <- NULL
  } else {   
    fill <- function(x, x_n, n_av, how = c("NA", "0")) {
      n_unav <- setdiff(x_n, n_av)
      n_unav_index <- match(n_unav, x_n)
      update_x <- x
      if (how == "NA") update_x[n_unav_index] <- as.numeric(NA)
      if (how == "0") update_x[n_unav_index] <- 0
      return(update_x)
    }  
    if(sequential == TRUE){
      if(to_sample1[1] == "none") {
        avail_sp <- avail_sp
      } else{
        avail_sp[to_sample1] <- 1        
      }
      com_av <- com[ , as.logical(avail_sp), drop = FALSE]
    } 
    to_sample2 <- alply(as.matrix(relab_com), 1, sampler, avail_sp, thresh)
    attributes(to_sample2) <- NULL
    if(all(unlist(lapply(to_sample2, is.null)))) {
      out$Scenario2 <- NULL
      warning("Threshold reached by Scenario 1 alone", call. = FALSE)
    } else {
      com_av_spl <- split(com_av, 1:nrow(com_av))
      #com_av_spl <- lapply(com_av_spl, function(x) x <- x[!is.na(x)])
      com_av_spl_names <- lapply(com_av_spl, function(x) names(x) <- colnames(com_av))
      av_names <- mapply(c, com_av_spl_names, to_sample2, SIMPLIFY = FALSE)
      com_av2 <- t(mapply(fill, split(com, 1:nrow(com)), list(colnames(com)), 
                          av_names, MoreArgs = list(how = "NA")))
      storage.mode(com_av2) <- "numeric"
      dimnames(com_av2) <- dimnames(com)
      relab_com_av2 <- t(mapply(fill, split(relab_com, 1:nrow(com)), 
                                list(colnames(com)), av_names, MoreArgs = list(how = "0")))
      dimnames(relab_com_av2) <- dimnames(com)
      storage.mode(relab_com_av2) <- "numeric"
      tot_pool_av2 <- sum(com_av2, na.rm = TRUE)/sum(com)
      tot_com_av2 <- rowSums(relab_com_av2)
      to_sample2 <- lapply(to_sample2, function(x) {
        if(is.null(x)) {
          x = "none"
        } else{
          x = x
        }
      })
      out$Scenario2$sample_species <- lapply(to_sample2, sort)
      out$Scenario2$com_available <- com_av2
      out$Scenario2$available_pool <- tot_pool_av2
      out$Scenario2$available_plot <- tot_com_av2      
    }
  }
  return(out)
}
