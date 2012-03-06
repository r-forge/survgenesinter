\name{generatePatientData}
\alias{generatePatientData}
\title{Generate Patient Data}
\usage{
  generatePatientData(N, l.1.tick, lambda)
}
\arguments{
  \item{N}{the sample size, i.e. the number of patients}

  \item{l.1.tick}{the (anticipated) lenth of the
  recruitment period}

  \item{lambda}{the mean survival time of the patients}
}
\value{
  \item{arrivalTimes}{the vector of simulated arrival
  times} \item{survivalTimes}{the vector of simulated
  survival times} \item{N.1}{number of patients with
  survival time less then the mean survival time 'lambda'}
  \item{N.2}{number of patients with survival time greater
  then the mean survival time 'lambda'} \item{l.1.tick}{the
  length of the recruitment period (unchanged parameter)}
  \item{lambda}{the mean survival time (parameter)}
}
\description{
  This function simulates entry dates of the patient data
  (e.g. date of diagnosis od data of surgery) and survival
  times of the patients.
}
\examples{
## the number of patients
N <- 50

## the anticipated length of the recruitment period
l.1.tick <- 60

## the mean survival time
lambda <- 60

\dontrun{P <- generatePatientData(N, l.1.tick, lambda)}
}
\author{
  Andreas Leha \email{andreas.leha@med.uni-goettingen.de}
}
\seealso{
  \code{\link{generateExpressionData}}
}

