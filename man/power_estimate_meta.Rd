\name{power_estimate_meta}
\alias{power_estimate_meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Calculation for Meta-analysis
}
\description{
This function calculates power for given parameters
}
\usage{
power_estimate_meta(k, es, total_n, heterogeneity = "fixed", alpha = .05, two_tailed = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{The number of studies included in the analysis}
  \item{es}{The minimum effect size of interest, expressed as the standardized mean difference}
  \item{total_n}{the total number of individuals in each comparison (both control and experimental groups)}
  \item{heterogeneity}{The level of heterogeneity expected. Can take the values of "small", "moderate", "large", or "fixed." If the meta-analysis is a fixed-effects meta-analysis, use "fixed" for heterogeneity.}
  \item{alpha}{The significance level of the meta-analysis}
  \item{two_tailed}{Whether the analysis will conduct a two-tailed test}
}
\details{
A power analysis is used to estimate the probability of finding a true effect of a given size. Researchers can conduct a power analysis after performing a scoping review to determine whether there are enough studies to perform an analysis. They can also conduct a power analysis after conducting a meta-analysis to determine whether they would have found an effect of interest with their number of studies. This function uses formulas provided by Valentine et al (2010) to estimate power for a given scenario.
}
\value{
This function returns the expected power probability
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
# Example 1: A one-tailed fixed-effects model including 40 studies, with 20 total individuals per comparison, and an anticipated effect size of .15
power_estimate_meta(k = 40, es = .15, total_n = 20, heterogeneity = "fixed", alpha = .05, two_tailed = FALSE)

# Example 2: A one-tailed random-effects model including 40 studies, with 40 total individuals per comparison, an anticipated effect size of .15, and moderate heterogeneity
power_estimate_meta(k = 40, es = .15, total_n = 20, heterogeneity = "fixed", alpha = .05, two_tailed = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
