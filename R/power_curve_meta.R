#' Power Curve Generation for Meta-analysis
#'
#' This function generates a power curve for a meta-analysis
#'
#' @usage power_curve_meta(lower_limit = 1, upper_limit = 100, es, total_n, heterogeneity = "fixed", alpha =.05, two_tailed = TRUE)
#'
#' @param lower_limit The lower bound of the plot
#' @param upper_limit The upper bound of the plot
#' @param es The minimum effect size of interest, expressed as the standardized mean difference
#' @param total_n the total number of individuals in each comparison (including control and experimetnal groups)
#' @param heterogeneity the level of heterogeneity expected. Can take the values of "small", "moderate", "large",
#' or "fixed." If the meta-analysis is a fixed-effects meta-analysis, use "fixed" for heterogeneity.
#' @param alpha The significance level of the meta-analysis
#' @param two_tailed Whether the analysis will conduct a two-tailed test
#'
#' @details A power curve is used to determine the number of studies needed to find a true effect of a given size.
#' Researchers can use a power curve to determine whether there are enough studies to perform an analysis.
#' They can also use a power curve after conducting a meta-analysis to determine whether they would have found an
#' effect of interest with their number of studies. This function uses formulas provided by Valentine et al (2010)
#' to estimate power for a given scenario.
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
#' # Example 1: A power curve for a one-tailed fixed-effects model. We assume we will have 20 total individuals
#' per comparison, and we have an anticipated effect size of .15
#' power_curve_meta(es = .15, total_n = 20, heterogeneity = "fixed", alpha = .05, two_tailed = FALSE)
#'
#'
#'

power_curve_meta = function(lower_limit = 1,
                            upper_limit = 100,
                            es,
                            total_n,
                            heterogeneity_vector = "fixed",
                            alpha =.05,
                            two_tailed = TRUE,
                            col_list = c("black", "blue4", "blue1", "turquoise")
                            ){

  if (missing(es)) {
    stop("The effect size of interest (es) must be included")
  }
  if(es < 0) {
    stop("The effect size (es) must be positive")
  }
  if(total_n < 0) {
    stop("The number of participants must be positive")
  }
  if(alpha < 0 | alpha > 1){
    stop("The confidence level (alpha) must be between 0 and 1")
  }
  if(lower_limit < 0) {
    stop("The lower bound of the graph must be greater than zero")
  }
  if(upper_limit <= lower_limit) {
    stop("The upper bound of the graph must be greater than the lower limit")
  }
  if(missing(total_n)) {
    stop("The total number of participants per comparison (total_n) must be included")
  }

  heterogeneity <- heterogeneity_vector[1]

  if (heterogeneity == "small") {
    hg_num <- .33
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

  power_vector <- numeric(length = upper_limit - lower_limit + 1)

  for (i in 1:length(power_vector)) {

    k <- lower_limit + i - 1

    var_individual_es <- total_n / ((total_n / 2)^2) + (es^2) / (2 * (total_n))
    rand_effect_correction <- hg_num * var_individual_es #if no heterogeneity, this will be 0
    corrected_var_ind_es <- rand_effect_correction + var_individual_es
    var_overall_es <- corrected_var_ind_es / k
    lamda <- es/sqrt(var_overall_es)
    power_vector[i] <- (1 - pnorm(z_needed - lamda))
  }



  plot(lower_limit:upper_limit, power_vector,
       xlab = "Number of Studies",
       ylab = "Power",
       main = "Power Curve",
       sub = paste("assuming ES of ", es, ", N of ", total_n, ", and ", paste(heterogeneity_vector, collapse = ", "), " heterogeneity", sep = ""),
       type = "l")


  if (length(heterogeneity_vector) > 1){
    for (het_val in 2:length(heterogeneity_vector)){
      heterogeneity <- heterogeneity_vector[het_val]

      if (heterogeneity == "small") {
        hg_num <- .33
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

      for (i in 1:length(power_vector)) {

        k <- lower_limit + i - 1

        var_individual_es <- total_n / ((total_n / 2)^2) + (es^2) / (2 * (total_n))
        rand_effect_correction <- hg_num * var_individual_es #if no heterogeneity, this will be 0
        corrected_var_ind_es <- rand_effect_correction + var_individual_es
        var_overall_es <- corrected_var_ind_es / k
        lamda <- es/sqrt(var_overall_es)
        power_vector[i] <- (1 - pnorm(z_needed - lamda))
      }

      lines(x = lower_limit:upper_limit, y = power_vector, col = col_list[het_val])

    }
    legend(x = "bottomright",
           legend = heterogeneity_vector,
           fill = col_list[1: length(heterogeneity_vector)],
           title = "Heterogeneity")
  }

}
