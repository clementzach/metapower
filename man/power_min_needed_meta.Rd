\name{power_min_needed_meta}
\alias{power_min_needed_meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum number of studies needed for Meta-analysis
}
\description{
This function calculates the minimum number of studies to generate a specificed level of power
}
\usage{
power_min_needed_meta(es, total_n, heterogeneity = "fixed", alpha = .05, two_tailed = TRUE, power_level = .80)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{The minimum effect size of interest, expressed as the standardized mean difference}
  \item{total_n}{the total number of individuals in each comparison (both control and experimental groups)}
  \item{heterogeneity}{The level of heterogeneity expected. Can take the values of "small", "moderate", "large", or "fixed." If the meta-analysis is a fixed-effects meta-analysis, use "fixed" for heterogeneity.}
  \item{alpha}{The significance level of the meta-analysis}
  \item{two_tailed}{Whether the analysis will conduct a two-tailed test}
  \item{power_level}{The desired power level. This is generally .80}
}
\details{This function calculates the minimum number of studies required to obtain a desired level of power. Researchers can conduct a power analysis after performing a scoping review to determine whether there are enough studies to perform an analysis. This function uses formulas provided by Valentine et al (2010) to determine the minimum number of studies required to achieve a desired level of power in a given scenario. A power level of .80 is conventionally desired for experiments, but in meta-analyses, it is more difficult to change the number of studies included in an analysis.

}
\value{
An integer which is the minimum number of studies for power to be above the specificed level.
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
# Example 1: The number of studies needed to achieve 80 percent power using a one-tailed fixed-effects model.
We assume we will have 20 total individuals per comparison, and we have an anticipated effect size of .15
power_min_needed_meta(es = .15, total_n = 20, heterogeneity = "fixed", two_tailed = FALSE, power_level = .80)

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
