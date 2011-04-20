\name{createMeanResultTable}
\alias{createMeanResultTable}
\title{Restructure the Simulation Output for Display...}
\usage{createMeanResultTable(resultTable, M.1, M.2, numSimRuns=1)}
\description{Restructure the Simulation Output for Display}
\details{This function restructures the output of \code{\link{interimTrial}} when run in
standalone mode and is also called internally by
\code{\link{interimTrialSimulation}}.}
\value{the parameter resultTable restructured as matrix}
\author{Andreas Leha \email{andreas.leha@med.uni-goettingen.de}}
\arguments{\item{resultTable}{the vector to bring in matrix form}
\item{M.1}{the number of analyses during the recruitment period to be simulated}
\item{M.2}{the number of analyses during the follow-up period to be simulated}
\item{numSimRuns}{the number of survival studies to be simulated}}
