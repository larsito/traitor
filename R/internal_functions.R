#' Fix total relative abundance to a set value
#' 
#' Takes a single plot (sample) and omits species, starting with the least 
#' abundant, until a user-set total relative abundance in \code{plot} is reached
#' @param plot a vector of species abundances
#' @param remain the total relative abundance value that \code{plot} should be 
#'   reduced to
#' @return a list of two vectors of species abundances, the first with the 
#'   original data provided by \code{plot}, the second where the total relative 
#'   abundance is that set by \code{remain}
#' @keywords internal

fixAbundance <- function(plot, remain) {
  p <- plot
  red <- 1 - remain
  l = vector("list", length = 2)
  names(l) <- c("Original", "Reduced")
  if (!is.vector(p)) stop("sample needs to be a vector")
  p <- p[p > 0]
  p <- p/sum(p)
  p <- sort(p)
  p_cs <- cumsum(p)
  p_u <- unique(p)
  l[[1]] <- p
  ll <- length(l)
  
# randomize the order of species that have exactly the same proportinal
# abundance, which would otherwise be alphabetical
# !! this leads to different species being excluded when applied with a 
# vector of remainAbu. So maybe better omit this (or make optional).
#   if (length(p_u) != length(p)) {
#     ma <- match(p, p_u)
#     name <- vector("character", length(p))
#     for(j in 1:length(unique(ma))){      
#       samp <- sample(names(p[which(ma==j)]))
#       name[which(ma==j)] <- samp
#     }
#     names(p) <- name
#     p_cs <- cumsum(p)
#   }
   
  if (sum(p) == remain) l[2] <- p
  if (sum(p) > remain) {
    if (any(p_cs == red)) {
      remov <- 1:which(p_cs == red)[1]
      p <- p[-remov]
      l[[2]] <- p    
    } else{
      if (all(p_cs > red)) { 
        p[1] <- p[1] - red
        l[[2]] <- p
      } else{
        remov <- 1:length(which(p_cs < red))
        p_r <- p_cs[length(remov)]
        p <- p[-remov]
        p[1] <- p[1] - (red - p_r)  
        l[[2]] <- p
      }    
    } 
  } else {
    stop("remain is set higher than the summed relative abundance in plot")
  }
return(l)
}

#' reduce abundance of species in community matrix
#' 
#' reduces the abundance of species in a community matrix according to a sample and a pool wise scenario for reduction
#' @param com a community data matrix with species as columnes and plots (samples) as rows.
#' @param remain the total relative abundance value that com should be reduced to.
#' @return a list of community data matrices: (1) Original community matrix, (2) Community matrix reduced according to Scenario 1, (3) Community matrix reduced according to Scenario 2
#' @keywords internal
reducer <- function(com, remain = 0.5) {
  if(is.matrix(com)) com <- as.data.frame(com)
  L <- list() # list to store results (i.e. different reduced community matrices)
  # reduce species on the basis of overall site relative abundances
  rel.ab.pool <- colSums(com)/sum(com)
  rel.ab.pool <- sort(rel.ab.pool)
  rel.ab.com <- com/rowSums(com)
  
  ## Scenario 1
  lr <- fixAbundance(rel.ab.pool, remain)  
  keep <- names(lr[[2]])
  out <- setdiff(names(rel.ab.pool), keep)
  com.red <- rel.ab.com[, -(na.omit(match(out, names(com))))]
  #fill <- rowSums(com.red)
  
  #   ## Scenario 2a
  #   lrl <- apply(com.red, 1, fixAbundance, remainAbu = min(fill), rel = F)
  #   lrl.low <- lapply(lrl, "[[", 2) 
  #   com.thresh.low <- matrix(0, nrow(com), ncol(com), dimnames = 
  #                                dimnames(com))
  #   
  #   for (i in 1:nrow(com)){  
  #     name  <- colnames(com.thresh.low)
  #     r <- lrl.low[[i]]
  #     index <- na.omit(match(names(r), name))
  #     com.thresh.low[i , index] <- r
  #   }
  
  # Scenario 2
  lrf <- apply(com, 1, fixAbundance, remain = remain)
  lrf.fix <- lapply(lrf, "[[", 2) 
  com.thresh.fix <- matrix(0, nrow(com), ncol(com), dimnames = 
                             dimnames(com))
  
  for (i in 1:nrow(com)){  
    name  <- colnames(com.thresh.fix)
    r <- lrf.fix[[i]]
    index <- na.omit(match(names(r), name))
    com.thresh.fix[i , index] <- r
  }
  
  L[[1]] <- rel.ab.com
  L[[2]] <- com.red
  #L[[3]] <- out(com.thresh.low)
  L[[3]] <- as.data.frame(com.thresh.fix[rowSums(com.thresh.fix)!=0, colSums(com.thresh.fix)!=0, drop = F])
  names(L) <- c("Orig", "Sc1", "Sc2")
  return(L)
}

#' define species to be sampled
#' 
#' defines what species in a plot (sample) should be sampled for measurement,
#' given what species have available trait data and how much of the relative
#' plot abundance should be covered
#' @param plot a vector of species abundances
#' @param avail.sp a binary vector of 0 and 1 with the same number and order of
#'   species as in \code{com} where 0 means that no trait data is available for
#'   that species
#' @param thresh a threshold value defining the minimum amount of relative 
#'   abundance that trait data should be available for
#' @return A character vector with the names of species that need to be sampled
#' @keywords internal
#' @export
sampler <- function(plot, avail_sp, thresh) {
  p <- plot
  av <- avail_sp
  if(any(is.na(av))) av <- as.numeric(!is.na(av))
  p_s <- sort(p)  # sort species by relative abundance
  if (is.null(names(av))) {
    names(av) <- names(p)
    warning("avail_sp is without names. Species are assumed to be in the same order as in plot.")
  }
  av_s <- av[names(p_s)]  # sort av to be in the same order
  p_av <- p_s[av_s == 1]  # vector of available species
  p_unav <- p_s[av_s == 0]  # vector of unavailable species
  p_av_tot <- sum(p_av)  # total abundance of available species
  samp <- vector("character")
  if (p_av_tot >= thresh) {
    samp <- NULL
  } else {
    p_unav_cs <- cumsum(rev(p_unav))
    add <- p_av_tot + p_unav_cs
    suff <- which(add > thresh)[1]
    samp <- names(rev(p_unav))[1:suff]
  }
  return(samp)
} 


#' wrapper for FD::dbFD 
#'
#' subsets the community and trait data to contain the same species, as needed by dbFD
#' @param com a community data matrix with species as columnes and plots (samples) as rows
#' @param trait vector, matrix, or data frame of trait values
#' @return a list of different indices of functional diversity. See FD::dbFD for
#'    details.
#' @import FD
#' @keywords internal
functDiv <- function(com, trait) {
  ins <- intersect(colnames(com), rownames(trait))
  com <- com[, ins, drop = FALSE]
  trait <- trait[ins, , drop = FALSE]
  out <- dbFD(trait, com, messages = FALSE, calc.FDiv = FALSE)
  return(out)
} 