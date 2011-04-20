\name{plotgendat}
\alias{plotgendat}
\title{This function visualizes the differential expression (fold changes) of the genes in...}
\usage{plotgendat(gendat)}
\description{This function visualizes the differential expression (fold changes) of the genes in
gendat among the groups of long survivers and short survivers.}
\details{The fucntion \code{generateExpressionData} has as one parameter the log fold change
between the gene expression levels of the group of long surviver and
the group of short survivers.  This parameter is visualized in this
function with a barplot.}
\value{none}
\author{Andreas Leha \email{andreas.leha@med.uni-goettingen.de}}
\seealso{\code{generateExpressionData}, \code{plotgeneratePatientData}}
\arguments{\item{gendat}{of type 'list'.  The return value of the function \code{generateExpressionData}}}
\examples{fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
Sigma.1 <- diag(1000)
Sigma.2 <- diag(1000)
N.1 <- 30
N.2 <- 30

gendat <- generateExpressionData(fc, Sigma.1, Sigma.2, N.1, N.2)
plotgendat(gendat)}
