\docType{package}
\name{survGenesInterim-package}
\alias{survgenesinterim-package}
\alias{survGenesInterim-package}
\title{Simulation of Survival Studies with Microarray Data}
\description{
  The main part of this package is the simulation of a
  survival study based on simulated gene expression level
  data and simulated patient data. Functions to generate
  such gene expression level data and patient data are
  included as well.  Resulting error rate (FDR) and power
  rate can be visualized. Thus, this package can be used
  for sample size or power estimations in the planning of a
  study.
}
\details{
  \tabular{ll}{ Package: \tab survGenesInterim\cr Type:
  \tab Package\cr Version: \tab 1.0\cr Date: \tab
  2010-07-19\cr License: \tab GPL (>= 2)\cr LazyLoad: \tab
  yes\cr }

  A frequent objective of clinical studies is to detect
  genes or biomarkers that can predict the outcome of
  therapy and thus the survival of patients.  Therefore,
  gene expression levels in samples of pretherapeutic ill
  or healthy samples are measured by DNA microarrays and
  compared to survival data of the patients.  Because
  usually, tissue samples are only available at
  distinguished time points and it can take several years
  until the full planned sample size is available and
  follow-up is complete.  In such long-lasting studies it
  is beneficial to obtain interim results already before
  their planned end. The details are of the simulation are
  given in Leha et al. (2010).

  This package allows to simulate such situation and to
  visualize the results, which is useful during the
  planning of this kind of studies, i.e. to assist in
  sample size planning.
}
\examples{
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
                                          numSimRuns)}
}
\author{
  Andreas Leha \email{andreas.leha@med.uni-goettingen.de}
}
\references{
  Leha, Andreas and Beissbarth, Tim and Jung, Klaus (2010):
  "Sequential Interim Analyses of Survival Data in
  Microarray Experiments", submitted
}
\keyword{package}

