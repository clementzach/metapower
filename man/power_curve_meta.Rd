\name{power_curve_meta}
\alias{power_curve_meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Curve Generation for Meta-analysis
}
\description{
This function generates a power curve for a meta-analysis
}
\usage{
power_curve_meta(lower_limit = 1, upper_limit = 100, es, total_n, heterogeneity = "fixed", alpha =.05, two_tailed = TRUE, col_list = c("black", "blue4", "blue1", "turquoise"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lower_limit}{The lower bound of the plot}
  \item{upper_limit}{The upper bound of the plot}
  \item{es}{The minimum effect size of interest, expressed as the standardized mean difference}
  \item{total_n}{the total number of individuals in each comparison (both control and experimental groups)}
  \item{heterogeneity_vector}{The level of heterogeneity expected. Can take the values of "small", "moderate", "large", or "fixed." If the meta-analysis is a fixed-effects meta-analysis, use "fixed" for heterogeneity. A vector can be passed with multiple heterogeneity levels.}
  \item{alpha}{The significance level of the meta-analysis}
  \item{two_tailed}{Whether the analysis will conduct a two-tailed test}
  \item{col_list}{Which colors will be used if doing multiple curves on the same plot. Cannot be shorter than heterogeneity_vector.}
}
\details{
A power curve is used to determine the number of studies needed to find a true effect of a given size. Researchers can use a power curve to determine whether there are enough studies to perform an analysis. They can also use a power curve after conducting a meta-analysis to determine whether they would have found an effect of interest with their number of studies. This function uses formulas provided by Valentine et al (2010) to estimate power for a given scenario.
}
\value{
This returns a plot with the power curve
}
\references{
Valentine, J. C., Pigott, T. D., & Rothstein, H. R. (2010). How many studies do you need? A primer on statistical power for meta-analysis. \emph{Journal of Educational and Behavioral Statistics, 35}(2), 215-247.
}
\author{
Zachary Clement
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Example 1: A power curve for a one-tailed model, including curves for both fixed and small heterogeneity. We assume we will have 20 total individuals per comparison, and we have an anticipated effect size of .15
power_curve_meta(es = .15, total_n = 20, heterogeneity = c("fixed", "small"), alpha = .05, two_tailed = FALSE, col_list = c("firebrick", "blue"))

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
