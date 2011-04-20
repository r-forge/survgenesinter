\name{interimTrialSimulation}
\alias{interimTrialSimulation}
\title{Simulation of Several Survival Studies With the Same Parameters...}
\usage{interimTrialSimulation(fc, Sigma.1, Sigma.2, N, l.1.tick, l.2, lambda,
    M.1, M.2, alpha, powerThreshold, adjustment, numSimRuns,
    parallel=TRUE)}
\description{Simulation of Several Survival Studies With the Same Parameters}
\details{\code{interimTrialSimulation} performs several simulations of survival
studies based on the given parameters.  For each survial study
simulation the patient data (arrival times and survival times) and the
gene expression level data are newly generated.}
\value{\item{resultTable}{the only real result. This table holds the mean error
and power rates from all the simulation runs.  The first row
additionally shows for each (interim) analysis the fraction of
simulated survival studies that were stopped then.}
The other values are copied versions of the corresponding input
parameters.  They are included into the result for convenience only.}
\author{Andreas Leha \email{andreas.leha@med.uni-goettingen.de}}
\seealso{calls \code{\link{generatePatientData}} and \code{\link{generateExpressionData}} to generate the
data and \code{interimTrial} to simulate a single survival study

\code{\link{plotinterimFDR}}, \code{\link{plotinterimAPR}}, and
\code{\link{plotinterimStops}} to visualize the results}
\arguments{\item{fc}{the vector of foldchanges between the two groups}
\item{Sigma.1}{the covariance matrix describing the correlation between the genes
in group one}
\item{Sigma.2}{the covariance matrix describing the correlation between the genes
in group two.  If this is NULL, the homoscedastic case is assumed
and Sigma.2 is set to Sigma.1.}
\item{N}{the sample size, i.e. the number of patients}
\item{l.1.tick}{the (anticipated) lenth of the recruitment period}
\item{l.2}{the length of the follow-up period to be simulated}
\item{lambda}{the mean survival time of the patients}
\item{M.1}{the number of analyses during the recruitment period to be simulated}
\item{M.2}{the number of analyses during the follow-up period to be simulated}
\item{alpha}{the error level at which the FDR is to be controlled}
\item{powerThreshold}{the study is stopped as soon as the estimated power rate
exceeds this threshold}
\item{adjustment}{the method to use for the p-value adjustment to account
for multiple testing.
In c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")}
\item{numSimRuns}{the number of survival studies to be simulated}
\item{parallel}{boolean value specifying whether to use the package 'multicore' for
a parallel execution of the simulation loop}}
\examples{## parameters controlling the gene level expression
fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
Sigma.1 <- diag(1000)
Sigma.2 <- diag(1000)

## parameters controlling the patient data
N <- 50
l.1.tick <- 60
l.2 <- 60
lambda <- 60

## parameters controlling the study design
M.1 <- 2
M.2 <- 2
alpha <- 0.05
powerThreshold <- 0.8
adjustment <- "BH"

## the number of studies to simulate
numSimRuns <- 2

result <- interimTrialSimulation(fc, Sigma.1, Sigma.2,
N, l.1.tick, l.2, lambda,
M.1, M.2, alpha, powerThreshold, adjustment,
numSimRuns, parallel=FALSE)
result$resultTable}
