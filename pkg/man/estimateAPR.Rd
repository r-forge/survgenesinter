\name{estimateAPR}
\alias{estimateAPR}
\title{Simple Estimation of the Average Power Rate}
\usage{
  estimateAPR(alpha, p_values, p_values_adjusted,
    lambda = 0.5)
}
\arguments{
  \item{alpha}{the type I error level}

  \item{p_values}{the raw p-values from the tests at level
  alpha}

  \item{p_values_adjusted}{the p-values adjusted for
  multiple testing}

  \item{lambda}{a tuning parameter.  Default value 0.5 as
  suggested by Storey and Tibshirany (2003)}
}
\value{
  The estimated APR, which is a single value in [0,1]
}
\description{
  Function for estimating the Average Power Rate (APR),
  also called Sensitivity.  Similar to the FDR estimator of
  Storey and Tibshirani (2003).
}
\examples{
## generate some p-values
p_values <- runif(100)
p_values_adjusted <- p.adjust(p_values, method="BH")

## estimate the power
\dontrun{estimateAPR(alpha=0.05, p_values=p_values, p_values_adjusted=p_values_adjusted)}
}
\author{
  Andreas Leha \email{andreas.leha@med.uni-goettingen.de}
}
\references{
  Storey, J. D. and R. Tibshirani (2003): "Statistical
  significance for genomewide studies", Proceedings of the
  National Academy of Sciences of the United States of
  America, 100, 9440-9445.
}
\seealso{
  \code{\link{estimateFDR}}
}

