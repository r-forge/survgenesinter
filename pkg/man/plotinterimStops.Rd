\name{plotinterimStops}
\alias{plotinterimStops}
\title{This function creates a barplot showing at which interim analyses the
simulated survival studies were stopped.}
\usage{
  plotinterimStops(interimdat)
}
\arguments{
  \item{interimdat}{of type 'list'.  The return value of
  either \code{interimTrial} or
  \code{interimTrialSimulation}}
}
\value{
  none
}
\description{
  For each simulated interim analysis one bar is drawn
  that's height represents the fraction of simulated
  survival studies, that were stopped at that particular
  interim analysis.
}
\examples{
## generate the data
fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
Sigma.1 <- diag(1000)
Sigma.2 <- diag(1000)

N <- 50
l.1.tick <- 60
l.2 <- 60
lambda <- 60

M.1 <- 2
M.2 <- 2
alpha <- 0.05
powerThreshold <- 0.8
adjustment <- "BH"

numSimRuns <- 2

\dontrun{result <- interimTrialSimulation(fc, Sigma.1, Sigma.2,
                                          N, l.1.tick, l.2, lambda,
                                          M.1, M.2, alpha, powerThreshold, adjustment,
                                          numSimRuns, parallel=FALSE)}

\dontrun{plotinterimStops(result)}
}
\author{
  Andreas Leha \email{andreas.leha@med.uni-goettingen.de}
}
\seealso{
  \code{plotinterimFDR}, \code{plotinterimAPR} for other
  plotting functions on the same data

  \code{interimTrial} and \code{interimTrialSimulation} for
  the data generation
}

