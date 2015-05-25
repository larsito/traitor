#' traitor: Tools For Functional Diversity Assessment With Missing Trait Data
#' 
#' The traitor package provides two main functionalities.
#' @section Evaluation of trait data availability for community data: In studies
#'   of trait based community ecology, researchers often face the problem that 
#'   they have a set of community data but do not have trait data for all 
#'   species. The \code{link{sampleSpecies}} function helps to assess how much 
#'   trait information is available for the community data and suggests which 
#'   species should be sampled to attain a more complete set of trait data.
#' @section Evaluation of sensitivity of FD indices to missing trait data: 
#'   \code{link{sampleSpecies}} artificially removes species from the community 
#'   data to demonstrate the effect of species with missing trait data on 
#'   measures of functional diversity. Unlike \code{link{sampleSpecies}}, this
#'   is not so much a practical tool for sampling campaigns, but providing a
#'   tool to investigate the issue of missing trait data on theoretical grounds.
#'   Generally, the function allows to analyse individual data sets as in
#'   Pakeman and Quested (2007), Pakeman (2014), and Majekova et al.
#'   (unpublished). With the \code{plot_reduceFD} function the results of such
#'   analyses can be visualized. The generated plots show the decrease in
#'   correlation between reduced and original functional diversity indices with
#'   increasing amounts of species removed from the community.
#'   
#' @references Pakeman, R. J. and Quested, H. M. (2007) Sampling plant 
#'   functional traits: What proportion of the species need to be measured? 
#'   Applied Vegetation Science 10, 91-96.
#'   
#'   Pakeman, R. J. (2013) Functional trait metrics are sensitive to the 
#'   completeness of the species' trait data? Methods in Ecology and Evolution 
#'   5, 9-15.
#'   
#'   Majekova, M., et al. Evaluating functional diversity: missing trait data and
#'   the importance of species abundance structure and data transformation. 
#'   Unpublished.
#'   
#' @docType package
#' @name traitor
NULL
