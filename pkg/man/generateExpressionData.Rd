\name{generateExpressionData}
\alias{generateExpressionData}
\title{Generate Gene Expression Data...}
\usage{generateExpressionData(fc=rep(0, 100), Sigma.1=diag(100), Sigma.2,
    N.1=10, N.2=10, use_cholesky=FALSE)}
\description{Generate Gene Expression Data}
\details{This function simulates gene expression data based on the multivariate
normal distribution for two groups of samples.}
\value{\item{X.1}{the simulated gene expression levels of group one}
\item{X.2}{the simulated gene expression levels of group two}
\item{d}{the dimension, i.e. the number of genes}
\item{fc}{the fold change vector.  This is the unchanged parameter to
the function.}}
\author{Andreas Leha \email{andreas.leha@med.uni-goettingen.de}}
\seealso{\code{\link{generatePatientData}}}
\arguments{\item{fc}{the vector of foldchanges between the two groups}
\item{Sigma.1}{the covariance matrix describing the correlation between the genes
in group one}
\item{Sigma.2}{the covariance matrix describing the correlation between the genes
in group two.  If this is NULL, the case of equal covariances is assumed
and Sigma.2 is set to Sigma.1.}
\item{N.1}{the sample size of group one}
\item{N.2}{the sample size of group two}
\item{use_cholesky}{this is a boolean parameter that indicates whether the
covariance matrices are cholesky decomposed.  This is an enourmous speed up
when simulating.}}
\examples{## create a vector of fold changes
fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))

## uncorrelated genes
Sigma.1 <- diag(1000)
Sigma.2 <- diag(1000)

## the sample sizes
N.1 <- 30
N.2 <- 30

G <- generateExpressionData(fc, Sigma.1, Sigma.2, N.1, N.2)}
