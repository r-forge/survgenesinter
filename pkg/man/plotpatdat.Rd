\name{plotpatdat}
\alias{plotpatdat}
\title{This function visualizes the patient data that is generated with
\code{generatePatientData}.}
\usage{
  plotpatdat(patdat)
}
\arguments{
  \item{patdat}{of type 'list'.  The return value of the
  function \code{generatePatientData}}
}
\value{
  none
}
\description{
  The function \code{generatePatientData} generates arrival
  times and survival times of patients included in the
  simulated survival study.  This function visualizes this
  data with timeline segments.  The distinction in long
  survivers and short survivers is color coded.
}
\examples{
N <- 50
l.1.tick <- 60
lambda <- 60

patdat <- generatePatientData(N, l.1.tick, lambda)
plotpatdat(patdat)
}
\author{
  Andreas Leha \email{andreas.leha@med.uni-goettingen.de}
}
\seealso{
  \code{generatePatientData},
  \code{plotgenerateExpressionData}
}

