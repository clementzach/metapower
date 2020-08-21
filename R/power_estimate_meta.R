#' Power Calculation for Meta-analysis
#'
#' This function calculates power for given parameters
#'
#' @usage power_estimate_meta(k, es, total_n, heterogeneity = "fixed", alpha = .05, two_tailed = TRUE)
#'
#' @param k The number of studies included in the analysis
#' @param es The minimum effect size of interest, expressed as the standardized mean difference
#' @param total_n the total number of individuals in each comparison (including control and experimetnal groups)
#' @param heterogeneity the level of heterogeneity expected. Can take the values of "small", "moderate", "large",
#' or "fixed." If the meta-analysis is a fixed-effects meta-analysis, use "fixed" for heterogeneity.
#' @param alpha The significance level of the meta-analysis
#' @param two_tailed Whether the analysis will conduct a two-tailed test
#'
#' @details A power analysis is used to estimate the probability of finding a true effect of a given size.
#' Researchers can conduct a power analysis after performing a scoping review to determine whether there
#' are enough studies to perform an analysis. They can also conduct a power analysis after conducting a meta-analysis
#' to determine whether they would have found an effect of interest with their number of studies. This function uses formulas
#' provided by Valentine et al (2010) to estimate power for a given scenario.
#'
#'
#'
#' @references
#'
#'Valentine, J. C., Pigott, T. D., & Rothstein, H. R. (2010). How many studies do you need? A primer on statistical
#'power for meta-analysis. \emph{Journal of Educational and Behavioral Statistics, 35}(2), 215-247.
#'
#'
#' @author Zachary Clement
#'
#'
#'
#' @examples
#'
#' # Example 1: A one-tailed fixed-effects model including 40 studies, with 20 total individuals per comparison, and an
#' anticipated effect size of .15
#' power_estimate_meta(k = 40, es = .15, total_n = 20, heterogeneity = "fixed", alpha = .05, two_tailed = FALSE)
#'
#' # Example 2: A one-tailed random-effects model including 40 studies, with 40 total individuals per comparison, an
#' anticipated effect size of .15, and moderate heterogeneity
#' power_estimate_meta(k = 40, es = .15, total_n = 20, heterogeneity = "fixed", alpha = .05, two_tailed = FALSE)
#'
#'

power_estimate_meta = function(k,
                               es,
                               total_n,
                               heterogeneity = "fixed",
                               alpha = .05,
                               two_tailed = TRUE){
  if (missing(k)){
    stop("The number of studies (k) must be included")
  }
  if (missing(es)) {
    stop("The effect size of interest (es) must be included")
  }
  if(es < 0 | es > 1) {
    stop("The effect size (es) must be between zero and one")
  }
  if(missing(total_n)) {
    stop("The total number of participants per comparison (total_n) must be included")
  }

  if (heterogeneity == "small") {
    hg_num = .33
  }
  if (heterogeneity == "moderate") {
    hg_num <- 1.0
  }
  if (heterogeneity == "large") {
    hg_num <- 3.0
  }
  if (heterogeneity == "fixed") {
    hg_num <- 0
  }
  if (match(heterogeneity, c("small", "moderate", "large", "fixed"), nomatch = 0) == 0) {
    stop("Error: Heterogeneity string is not recognized. Use small, medium, large, or fixed")
  }

  if (two_tailed) {
    z_needed <- qnorm(1 - (alpha / 2))
  }
  else {
    z_needed <- qnorm(1 - alpha)
  }


  var_individual_es <- (total_n) / ((total_n / 2)^2) + (es^2) / (2 * (total_n))
  rand_effect_correction <- hg_num*(var_individual_es) #if no heterogeneity, this will be 0
  corrected_var_ind_es <- rand_effect_correction + var_individual_es
  var_overall_es <- corrected_var_ind_es / k
  lamda <- es/sqrt(var_overall_es)
  power <- 1 - pnorm(z_needed - lamda)

  return(power)


}
