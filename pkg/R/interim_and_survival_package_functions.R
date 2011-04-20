##' The main part of this package is the simulation of a survival study
##' based on simulated gene expression level data and simulated patient data.
##' Functions to generate such gene expression level data and patient data are
##' included as well.  Resulting error rate (FDR) and power rate can be visualized.
##' Thus, this package can be used for sample size or power estimations in the planning of a study.
##'
##' \tabular{ll}{
##' Package: \tab survGenesInterim\cr
##' Type: \tab Package\cr
##' Version: \tab 1.0\cr
##' Date: \tab 2010-07-19\cr
##' License: \tab GPL (>= 2)\cr
##' LazyLoad: \tab yes\cr
##' }
##'
##' A frequent objective of clinical studies is to detect genes or
##' biomarkers that can predict the outcome of therapy and thus the
##' survival of patients.  Therefore, gene expression levels in samples of
##' pretherapeutic ill or healthy samples are measured by DNA microarrays
##' and compared to survival data of the patients.  Because usually,
##' tissue samples are only available at distinguished time points and it
##' can take several years until the full planned sample size is available
##' and follow-up is complete.  In such long-lasting studies it is
##' beneficial to obtain interim results already before their planned end.
##' The details are of the simulation are given in Leha et al. (2010).
##'
##' This package allows to simulate such situation and to visualize the
##' results, which is useful during the planning of this kind of studies,
##' i.e. to assist in sample size planning.
##'
##' @name survGenesInterim-package
##' @aliases survgenesinterim-package
##' @docType package
##' @title Simulation of Survival Studies with Microarray Data
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @keywords package
##' @references Leha, Andreas and Beissbarth, Tim and Jung, Klaus (2010):
##'   "Sequential Interim Analyses of Survival Data in Microarray Experiments",
##'   submitted
##' @examples
##' fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
##' Sigma.1 <- diag(1000)
##' Sigma.2 <- diag(1000)
##' 
##' N <- 50
##' l.1.tick <- 60
##' l.2 <- 60
##' lambda <- 60
##' 
##' M.1 <- 2
##' M.2 <- 2
##' alpha <- 0.05
##' powerThreshold <- 0.8
##' adjustment <- "BH"
##' 
##' numSimRuns <- 2
##' 
##' result <- interimTrialSimulation(fc, Sigma.1, Sigma.2,
##'                                  N, l.1.tick, l.2, lambda,
##'                                  M.1, M.2, alpha, powerThreshold, adjustment,
##'                                  numSimRuns)
NA


##' Generate Patient Data
##'
##' This function simulates entry dates of the patient data (e.g. date of
##' diagnosis od data of surgery) and survival times of the patients.
##'
##' @param N the sample size, i.e. the number of patients
##' @param l.1.tick the (anticipated) lenth of the recruitment period
##' @param lambda the mean survival time of the patients
##' @return \item{arrivalTimes}{the vector of simulated arrival times}
##'  \item{survivalTimes}{the vector of simulated survival times}
##'  \item{N.1}{number of patients with survival time less then the mean
##'    survival time 'lambda'}
##'  \item{N.2}{number of patients with survival time greater then the mean
##'    survival time 'lambda'}
##'  \item{l.1.tick}{the length of the recruitment period (unchanged parameter)}
##'  \item{lambda}{the mean survival time (parameter)}
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @seealso \code{\link{generateExpressionData}}
##' @examples
##' ## the number of patients
##' N <- 50
##' 
##' ## the anticipated length of the recruitment period
##' l.1.tick <- 60
##' 
##' ## the mean survival time
##' lambda <- 60
##' 
##' P <- generatePatientData(N, l.1.tick, lambda)
generatePatientData <- function(N, l.1.tick, lambda)
{
  if(!is.numeric(N))
    stop("the sample size 'N' must be numeric")
  if(!is.vector(N))
    stop("the sample size 'N' must be a vector")
  if(!(length(N) == 1))
    stop("the sample size 'N' must be a single value")
  if(N <= 0)
    stop("the sample size 'N' must be positive")

  if(!is.numeric(l.1.tick))
    stop("the length of the recruitment period 'l.1.tick' must be numeric")
  if(!is.vector(l.1.tick))
    stop("the length of the recruitment period 'l.1.tick' must be a vector")
  if(!(length(l.1.tick) == 1))
    stop("the length of the recruitment period 'l.1.tick' must be a single value")
  if(l.1.tick <= 0)
    stop("the length of the recruitment period 'l.1.tick' must be positive")

  if(!is.numeric(lambda))
    stop("the mean survival time 'lambda' must be numeric")
  if(!is.vector(lambda))
    stop("the mean survival time 'lambda' must be a vector")
  if(!(length(lambda) == 1))
    stop("the mean survival time 'lambda' must be a single value")
  if(lambda <= 0)
    stop("the mean survival time 'lambda' must be positive")


  ## patients arrive uniformly distributed
  arrivalTimes <- runif(N, 0, l.1.tick)

  ## the survival times are exponetial distributed
  survivalTimes <- rexp(N, 1/lambda)

  N.1 <- sum(survivalTimes <= lambda)
  N.2 <- N - N.1
  
  patdat <- list(arrivalTimes=arrivalTimes,
                 survivalTimes=survivalTimes,
                 N.1=N.1,
                 N.2=N.2,
                 l.1.tick=l.1.tick,
                 lambda=lambda)
  return(patdat)
}


##' Generate Gene Expression Data
##'
##' This function simulates gene expression data based on the multivariate
##' normal distribution for two groups of samples.
##'
##' @param fc the vector of foldchanges between the two groups
##' @param Sigma.1 the covariance matrix describing the correlation between the genes
##'    in group one
##' @param Sigma.2 the covariance matrix describing the correlation between the genes
##'    in group two.  If this is NULL, the case of equal covariances is assumed
##'    and Sigma.2 is set to Sigma.1.
##' @param N.1 the sample size of group one
##' @param N.2 the sample size of group two
##' @param use_cholesky this is a boolean parameter that indicates whether the
##'    covariance matrices are cholesky decomposed.  This is an enourmous speed up
##'    when simulating.
##' @return   \item{X.1}{the simulated gene expression levels of group one}
##'   \item{X.2}{the simulated gene expression levels of group two}
##'   \item{d}{the dimension, i.e. the number of genes}
##'   \item{fc}{the fold change vector.  This is the unchanged parameter to
##'     the function.}
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @seealso \code{\link{generatePatientData}}
##' @examples
##' ## create a vector of fold changes
##' fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
##' 
##' ## uncorrelated genes
##' Sigma.1 <- diag(1000)
##' Sigma.2 <- diag(1000)
##' 
##' ## the sample sizes
##' N.1 <- 30
##' N.2 <- 30
##' 
##' G <- generateExpressionData(fc, Sigma.1, Sigma.2, N.1, N.2)
generateExpressionData <- function(fc=rep(0,100),
                   Sigma.1=diag(100),
                   Sigma.2=NULL,
                   N.1=10,
                   N.2=10,
                   use_cholesky=FALSE)
{
  if (!is.vector(fc))
    stop("the fold change vector 'fc' must be a vector")
  if(!is.numeric(fc))
    stop("the fold change vector 'fc' must be numeric")
  if(!is.matrix(Sigma.1))
    stop("the covariance matrix 'Sigma.1' must be a matrix")
  if(!is.numeric(Sigma.1))
    stop("the covariance matrix 'Sigma.1' must be numeric")
  if(!is.null(Sigma.2)) {
    if(!is.matrix(Sigma.2))
      stop("the covariance matrix 'Sigma.2' must be a matrix (or 'NULL')")
    if(!is.numeric(Sigma.2))
      stop("the covariance matrix 'Sigma.2' must be numeric (or 'NULL')")
  }
  if(!is.numeric(N.1))
    stop("the group size 'N.1' must be numeric")
  if(!is.vector(N.1))
    stop("the group size 'N.1' must be a vector")
  if(!(length(N.1) == 1))
    stop("the group size 'N.1' must be a single value")
  if(N.1 <= 0)
    stop("the group size 'N.1' must be positive")
  if(!is.numeric(N.2))
    stop("the group size 'N.2' must be numeric")
  if(!is.vector(N.2))
    stop("the group size 'N.2' must be a vector")
  if(!(length(N.2) == 1))
    stop("the group size 'N.2' must be a single value")
  if(N.2 <= 0)
    stop("the group size 'N.2' must be positive")
  if(!is.logical(use_cholesky))
    stop("the parameter 'use_cholesky' must be either 'TRUE' or 'FALSE'")

  if(length(fc) != dim(Sigma.1)[1])
    stop("dimensions of 'fc' and 'Sigma.1' do not match")
  if(!is.null(Sigma.2)) {
    if(!(length(dim(Sigma.1)) == length(dim(Sigma.2))))
      stop("dimensions of 'Sigma.1' and 'Sigma.2' do not match")
    if(!all((dim(Sigma.1) == dim(Sigma.2))))
      stop("dimensions of 'Sigma.1' and 'Sigma.2' do not match")
  }

  
  if(is.null(Sigma.2))
    Sigma.2 <- Sigma.1  
  
  d <- length(fc)

  mu.1 <- rep(0, d)
  mu.2 <- mu.1 + fc
  
  if (use_cholesky)
  {
    random_independent <- rnorm(N.1*d)
    dim(random_independent) <- c(N.1, d)
    X.1 <- random_independent %*% Sigma.1
    random_independent <- rnorm(N.2*d)
    dim(random_independent) <- c(N.2, d)
    X.2 <- random_independent %*% Sigma.2
  } else {
    X.1 <- rmvnorm(N.1, rep(0,d), Sigma.1)
    X.2 <- rmvnorm(N.2, rep(0,d), Sigma.2)
  }
  
  X.1 <- t(X.1)
  X.2 <- t(X.2)

  X.1 <- X.1 + mu.1
  X.2 <- X.2 + mu.2
  
  gendat <- list(X.1=X.1, X.2=X.2, d=d, fc=fc)
  return(gendat)
}


##' Simple Estimation of the Average Power Rate
##'
##' Function for estimating the Average Power Rate (APR), also called
##' Sensitivity.  Similar to the FDR estimator of Storey and Tibshirani (2003).
##'
##' @param alpha the type I error level
##' @param p_values the raw p-values from the tests at level alpha
##' @param p_values_adjusted the p-values adjusted for multiple testing
##' @param lambda a tuning parameter.  Default value 0.5 as suggested by
##'    Storey and Tibshirany (2003)
##' @return The estimated APR, which is a single value in [0,1]
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @references Storey, J. D. and R. Tibshirani (2003):
##'   "Statistical significance for genomewide studies",
##'   Proceedings of the National Academy of Sciences of the United States
##'   of America, 100, 9440-9445.
##' @seealso \code{\link{estimateFDR}}
##' @examples
##' ## generate some p-values
##' p_values <- runif(100)
##' p_values_adjusted <- p.adjust(p_values, method="BH")
##' 
##' ## estimate the power
##' estimateAPR(alpha=0.05, p_values=p_values, p_values_adjusted=p_values_adjusted)
estimateAPR <- function(alpha, p_values, p_values_adjusted, lambda=0.5)
{
  if(!is.numeric(alpha))
    stop("the error level 'alpha' must be numeric")
  if(!is.vector(alpha))
    stop("the error level 'alpha' must be a vector")
  if(!(length(alpha) == 1))
    stop("the error level 'alpha' must be a single value")
  if(alpha < 0 || alpha > 1)
    stop("the error level 'alpha' must be in [0,1]")

  if(!is.numeric(lambda))
    stop("the tuning parameter 'lambda' must be numeric")
  if(!is.vector(lambda))
    stop("the tuning parameter 'lambda' must be a vector")
  if(!(length(lambda) == 1))
    stop("the tuning parameter 'lambda' must be a single value")
  if(lambda < 0 || lambda > 1)
    stop("the tuning parameter 'lambda' must be in [0,1]")

  if (!is.vector(p_values))
    stop("'p_values' must be a vector")
  if(!is.numeric(p_values))
    stop("'p_values' must be numeric")

  if (!is.vector(p_values_adjusted))
    stop("'p_values_adjusted' must be a vector")
  if(!is.numeric(p_values_adjusted))
    stop("'p_values_adjusted' must be numeric")

  if (!(length(p_values) == length(p_values_adjusted)))
    stop("'p_values' and 'p_values_adjusted' must not differ in length")


  d <- length(p_values)

  pi_hat <- length(p_values[p_values > lambda]) / (d * (1-lambda))
  pi_hat <- min(pi_hat,1)

  p_values_ordered <- sort(p_values)
  alpha_adjusted <- 0
  for (i in 1:length(p_values_ordered))
  {
    if (p_values_ordered[i] <= (i/d)*alpha)
      alpha_adjusted <- (i/d)*alpha
  }
  
  d1_hat <- (1 - pi_hat) * d
  d0_hat <- pi_hat * d
  FP_hat <- alpha_adjusted * d0_hat
  R <- length(p_values_adjusted[p_values_adjusted<= alpha])
  TP_hat <- R - FP_hat

  if (d1_hat == 0) {
    0
  } else {
    TP_hat / d1_hat
  }
}


##' Simple Estimation of the False Discovery Rate
##'
##' Function for estimating the False Discovery Rate (FDR).
##' As in Storey and Tibshirani (2003).
##'
##' @param alpha the type I error level
##' @param p_values the raw p-values from the tests at level alpha
##' @param p_values_adjusted the p-values adjusted for multiple testing
##' @param lambda a tuning parameter.  Default value 0.5 as suggested by
##'    Storey and Tibshirany (2003)
##' @return The estimated FDR, which is a single value in [0,1]
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @references Storey, J. D. and R. Tibshirani (2003):
##'   "Statistical significance for genomewide studies",
##'   Proceedings of the National Academy of Sciences of the United States
##'   of America, 100, 9440-9445.
##' @seealso \code{\link{estimateAPR}}
##' @examples
##' ## generate some p-values
##' p_values <- runif(100)
##' p_values_adjusted <- p.adjust(p_values, method="BH")
##' 
##' ## estimate the FDR
##' estimateFDR(alpha=0.05, p_values=p_values, p_values_adjusted=p_values_adjusted)
estimateFDR <- function(alpha, p_values, p_values_adjusted, lambda=0.5)
{
  if(!is.numeric(alpha))
    stop("the error level 'alpha' must be numeric")
  if(!is.vector(alpha))
    stop("the error level 'alpha' must be a vector")
  if(!(length(alpha) == 1))
    stop("the error level 'alpha' must be a single value")
  if(alpha < 0 || alpha > 1)
    stop("the error level 'alpha' must be in [0,1]")

  if(!is.numeric(lambda))
    stop("the tuning parameter 'lambda' must be numeric")
  if(!is.vector(lambda))
    stop("the tuning parameter 'lambda' must be a vector")
  if(!(length(lambda) == 1))
    stop("the tuning parameter 'lambda' must be a single value")
  if(lambda < 0 || lambda > 1)
    stop("the tuning parameter 'lambda' must be in [0,1]")

  if (!is.vector(p_values))
    stop("'p_values' must be a vector")
  if(!is.numeric(p_values))
    stop("'p_values' must be numeric")

  if (!is.vector(p_values_adjusted))
    stop("'p_values_adjusted' must be a vector")
  if(!is.numeric(p_values_adjusted))
    stop("'p_values_adjusted' must be numeric")

  if (!(length(p_values) == length(p_values_adjusted)))
    stop("'p_values' and 'p_values_adjusted' must not differ in length")

  d <- length(p_values)
  
  pi_hat <- length(p_values[p_values > lambda]) / (d * (1-lambda)) ## might become > 1
  pi_hat <- min(pi_hat,1)

  p_values_ordered <- sort(p_values)
  alpha_adjusted <- 0
  for (i in 1:length(p_values_ordered))
  {
    if (p_values_ordered[i] <= (i/d)*alpha)
      alpha_adjusted <- (i/d)*alpha
  }
  
  d0_hat <- pi_hat * d
  FP_hat <- alpha_adjusted * d0_hat
  R <- length(p_values_adjusted[p_values_adjusted<= alpha])

  if (R == 0) {
    0
  } else {
    FP_hat / R
  }
}


##' \code{interimTrial} simulates one survival study with interim analyses based
##' on the patient data and the gene expression data it is supplied
##' with.  At each interim analyses the error rate (FDR) and the average power
##' rate (APR) are both calculated and estimated.
##'
##'   interimTrial comes in two variants: When \code{standalone} is 'TRUE'
##'   the result is nicely formatted and can be used with the plotting
##'   functions \code{plotinterimFDR}, \code{plotinterimAPR}, and
##'   \code{plotinterimStops}.  When \code{standalone} is 'FALSE' the result
##'   is presented as one single vector, which makes it usable within
##'   \code{interimTrialSimulation}.
##' 
##'   The \code{adjustment} parameter is used as parameter to
##'   \code{p.adjust}.
##' 
##'   The parameters \code{patdat} and \code{gendat} are of type \code{list}
##'   and do not necessarily have to be obtained from the functions
##'   \code{generatePatientData} and \code{generateExpressionData}.  Especially the entries \code{N.1}
##'   and \code{N.2} in \code{patdat} can be changed.  But in this case make
##'   sure to generate \code{gendat} accordingly.
##' 
##'   Internally interimTrial performs a gene-wise Cox-regression to
##'   determin the survival correlated genes.
##'
##' @param patdat the patient data (arrival times and survival times)
##'   as generated by \code{\link{generatePatientData}}
##' @param gendat the gene expression data as generated by \code{\link{generateExpressionData}}
##' @param l.2 the length of the follow-up period to be simulated 
##' @param M.1 the number of analyses during the recruitment period to be simulated
##' @param M.2 the number of analyses during the follow-up period to be simulated
##' @param alpha the error level at which the FDR is to be controlled
##' @param powerThreshold the study is stopped as soon as the estimated power rate
##'   exceeds this threshold
##' @param adjustment the method to use for the p-value adjustment to account
##'   for multiple testing.
##'   In c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
##' @param standalone boolean parameter that controls whether it is used within
##'   a bigger simulation.  Defaults to 'TRUE'
##' @return In case the function is used in standalone mode:
##'   \item{resultTable}{the only real result. This table holds the error
##'     and power rates that came from the simulation.  The fist row
##'   additionally shows at which analysis the simulated survival study was stopped.}
##'   The other values are copied versions of the corresponding input
##'   parameters.  They are included into the result for convenience only.
##' 
##'   In case the function is not used in standalone mode:
##'   The return value is one vector containing the elements of \code{resultTable}.
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @seealso \code{\link{generatePatientData}} and \code{\link{generateExpressionData}} to generate the input
##'   data
##' 
##'   \code{\link{plotinterimFDR}}, \code{\link{plotinterimAPR}}, and
##'   \code{\link{plotinterimStops}} to visualize the results
##' 
##'   \code{\link{interimTrialSimulation}} for a wrapper that simulates not
##'   one but many survival studies (and calls this function)
##' @examples
##' N <- 50
##' l.1.tick <- 60
##' lambda <- 60
##' 
##' patdat <- generatePatientData(N, l.1.tick, lambda)
##' 
##' fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
##' Sigma.1 <- diag(1000)
##' Sigma.2 <- diag(1000)
##' N.1 <- patdat$N.1
##' N.2 <- patdat$N.2
##' 
##' gendat <- generateExpressionData(fc, Sigma.1, Sigma.2, N.1, N.2)
##' 
##' l.2 <- 60
##' M.1 <- 2
##' M.2 <- 2
##' alpha <- 0.05
##' powerThreshold <- 0.8
##' adjustment <- "BH"
##' 
##' T <- interimTrial(patdat, gendat, l.2, M.1, M.2, alpha, powerThreshold, adjustment)
interimTrial <- function(patdat,
                         gendat,
                         l.2,
                         M.1,
                         M.2,
                         alpha,
                         powerThreshold,
                         adjustment,
                         standalone=TRUE)
{
  if(!is.list(patdat))
    stop("'patdat' must be the return value from the function 'generatePatientData'")
  if(!is.list(gendat))
    stop("'gendat' must be the return value from the function 'gendat'")

  if(!is.numeric(l.2))
    stop("the length of the follow-up period 'l.2' must be numeric")
  if(!is.vector(l.2))
    stop("the length of the follow-up period 'l.2' must be a vector")
  if(!(length(l.2) == 1))
    stop("the length of the follow-up period 'l.2' must be a single value")
  if(l.2 <= 0)
    stop("the length of the follow-up period 'l.2' must be positive")

  if(!is.numeric(M.1))
    stop("the number of analyses during recruitment period 'M.1' must be numeric")
  if(!is.vector(M.1))
    stop("the number of analyses during recruitment period 'M.1' must be a vector")
  if(!(length(M.1) == 1))
    stop("the number of analyses during recruitment period 'M.1' must be a single value")
  if(M.1 <= 0)
    stop("the number of analyses during recruitment period 'M.1' must be positive")
  
  if(!is.numeric(M.2))
    stop("the number of analyses during follow-up period 'M.2' must be numeric")
  if(!is.vector(M.2))
    stop("the number of analyses during follow-up period 'M.2' must be a vector")
  if(!(length(M.2) == 1))
    stop("the number of analyses during follow-up period 'M.2' must be a single value")
  if(M.2 <= 0)
    stop("the number of analyses during follow-up period 'M.2' must be positive")

  if(!is.numeric(alpha))
    stop("the error level 'alpha' must be numeric")
  if(!is.vector(alpha))
    stop("the error level 'alpha' must be a vector")
  if(!(length(alpha) == 1))
    stop("the error level 'alpha' must be a single value")
  if(alpha < 0 || alpha > 1)
    stop("the error level 'alpha' must be in [0,1]")

  if(!is.numeric(powerThreshold))
    stop("the 'powerThreshold' must be numeric")
  if(!is.vector(powerThreshold))
    stop("the 'powerThreshold' must be a vector")
  if(!(length(powerThreshold) == 1))
    stop("the 'powerThreshold' must be a single value")
  if(powerThreshold < 0 || powerThreshold > 1)
    stop("the 'powerThreshold' must be in [0,1]")

  if (!(adjustment %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                          "fdr", "none")))
    stop("the adjustment method 'adjustment' must be in c(\"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\", \"none\")")


  
  ## #####
  ## extract data
  ##
  arrivalTimes <- patdat$arrivalTimes
  survivalTimes <- patdat$survivalTimes
  lambda <- patdat$lambda
  d <- gendat$d
  ##
  N <- patdat$N.1 + patdat$N.2
  M <- M.1 + M.2
  ##
  real_num_differential <- sum(gendat$fc != 0)
  indexDifferentialGenes <- which(gendat$fc != 0)
  d.1 <- length(indexDifferentialGenes)
  d.0 <- d - d.1
  ##
  eventTimes <- arrivalTimes + survivalTimes
  ##
  survivalTimesSorted <- sort(patdat$survivalTimes)
  groupCutoff <- survivalTimesSorted[patdat$N.1]
  groupShort <- patdat$survivalTimes <= groupCutoff
  groupLong <- 1:N
  groupLong <- groupLong[!groupShort]
  ##  
  X <- numeric(d*N)
  dim(X) <- c(d, N)
  X[,groupShort] <- gendat$X.1
  X[,groupLong] <- gendat$X.2
  ## #####

  
  ## when in l.1 will interim analyses be conducted?
  arrivalTimesSorted <- sort(arrivalTimes)
  interimCutoffs <- c()
  for (m in 1:M.1)
  {
    interimCutoffs <-
      c(interimCutoffs,
        arrivalTimesSorted[length(arrivalTimesSorted)*(m/M.1)])
  }
  ## when in l.2 will interim analyses be conducted?
  analysisTimes <-
    c(interimCutoffs,
      (((1:M.2) * (l.2/M.2)) + interimCutoffs[length(interimCutoffs)]))

  
  ## which samples are available at each stage?
  interimIndex <- list()
  for (m in 1:M)
  {
    interimIndex[[m]] <- which(arrivalTimes <= analysisTimes[m])
  }
  
  resultVector <- c()
  stoppingAnalysis <- NA
  
  for (m in 1:M)
  {
    if (is.na(stoppingAnalysis)) {
      ## the survival times are only knwon until the time of the analysis
      survivalTimesCensored <- pmin(analysisTimes[m]-arrivalTimes,
                                    survivalTimes)

      
      ## calculation of the censor variable
      censor <- rep(0,N)
      censor[which(eventTimes < analysisTimes[m])] <- 1
      
      ## what data are avalable?
      survivalTimesCensored.currentAnalysis <- survivalTimesCensored[interimIndex[[m]]]
      censor.currentAnalysis <- censor[interimIndex[[m]]]
      X.currentAnalysis <- X[,interimIndex[[m]]]
      
      ## order the data
      currentAnalysis.order <- order(survivalTimesCensored.currentAnalysis)
      survivalTimesCensored.currentAnalysis <-
        round(survivalTimesCensored.currentAnalysis[currentAnalysis.order], 1)
      censor.currentAnalysis <- censor.currentAnalysis[currentAnalysis.order]
      X.currentAnalysis <- X.currentAnalysis[,currentAnalysis.order]
      
      ## p-value calculation by gene-wise cox regression
      P <- rep(NA, d)
      K <- Surv(survivalTimesCensored.currentAnalysis, censor.currentAnalysis)
      for (j in 1:d) {
        C <- try(summary(suppressWarnings(coxph(K ~ X.currentAnalysis[j,])))$logtest,
                 silent=TRUE)
        P[j] <- C[3]
      }
      P <- as.numeric(P)
      
      ## multiple adjustment
      P_adjusted <- p.adjust(P, method=adjustment)


      ## ########
      ## calulating rates
      ##
      ## the number of false positives
      if (real_num_differential == 0) {
        FP <- sum(P_adjusted[!is.na(P_adjusted)] < alpha)
      } else {
        FP <- sum(P_adjusted[-indexDifferentialGenes][!is.na(P_adjusted[-indexDifferentialGenes])] < alpha)
      }
      ##
      ## the number of rejected
      R <- sum(P_adjusted[!is.na(P_adjusted)]< alpha)
      ##
      ## the number of true positives
      TP <- R - FP
      ##
      ## ratio of false positives to recected = false positive proportion (FDP)
      if (R > 0) {
        FDP <- FP/R
      } else {
        FDP <- 0
      }
      ##
      ## ratio of false positives to recected = average power proportion (APP)
      if (d.1 > 0) {
        APP <- TP/d.1
      } else {
        APP <- 0
      }
      ## ########

      ## ########
      ## estimating rates
      ##
      ## estimated average power rate (eAPR)
      eAPP <- estimateAPR(alpha, P, P_adjusted)
      ##
      ## estimated false discovery rate (eFDR)
      eFDP <- estimateFDR(alpha, P, P_adjusted)
      ## ########


      ## early stopping ?
      if ((!is.na(eAPP)) && (eAPP > powerThreshold)) {
        stoppingAnalysis <- m
      }
      
      ## fill result vector
      resultVector <- c(resultVector,
                        FDP,
                        APP,
                        eFDP,
                        eAPP)

    } else {
      resultVector <- c(resultVector,
                        NA,
                        NA,
                        NA,
                        NA)
    }
  }
  if (is.na(stoppingAnalysis))
    stoppingAnalysis <- M
  resultVector <- c(stoppingAnalysis, resultVector)
  
  if (standalone) {
    resultTable <- createMeanResultTable(as.matrix(resultVector), M.1, M.2)
    return(list(fc=gendat$fc,
                lambda=patdat$lambda,
                N=N,
                l.1.tick=patdat$l.1.tick,
                l.2=l.2,
                M.1=M.1,
                M.2=M.2,
                alpha=alpha,
                powerThreshold=powerThreshold,
                adjustment=adjustment,
                numSimRuns=1,
                resultTable=resultTable))
  } else {
    return(resultVector)
  }
}


##' Simulation of Several Survival Studies With the Same Parameters
##'
##' \code{interimTrialSimulation} performs several simulations of survival
##' studies based on the given parameters.  For each survial study
##' simulation the patient data (arrival times and survival times) and the
##' gene expression level data are newly generated.
##'
##' @param fc the vector of foldchanges between the two groups
##' @param Sigma.1 the covariance matrix describing the correlation between the genes
##'   in group one
##' @param Sigma.2 the covariance matrix describing the correlation between the genes
##'   in group two.  If this is NULL, the homoscedastic case is assumed
##'   and Sigma.2 is set to Sigma.1.
##' @param N the sample size, i.e. the number of patients
##' @param l.1.tick the (anticipated) lenth of the recruitment period
##' @param l.2 the length of the follow-up period to be simulated
##' @param lambda the mean survival time of the patients
##' @param M.1 the number of analyses during the recruitment period to be simulated
##' @param M.2 the number of analyses during the follow-up period to be simulated
##' @param alpha the error level at which the FDR is to be controlled
##' @param powerThreshold the study is stopped as soon as the estimated power rate
##'   exceeds this threshold
##' @param adjustment the method to use for the p-value adjustment to account
##'   for multiple testing.
##'   In c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
##' @param numSimRuns the number of survival studies to be simulated
##' @param parallel boolean value specifying whether to use the package 'multicore' for
##'   a parallel execution of the simulation loop
##' @return   \item{resultTable}{the only real result. This table holds the mean error
##'     and power rates from all the simulation runs.  The first row
##'     additionally shows for each (interim) analysis the fraction of
##'     simulated survival studies that were stopped then.}
##'   The other values are copied versions of the corresponding input
##'   parameters.  They are included into the result for convenience only.
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @seealso   calls \code{\link{generatePatientData}} and \code{\link{generateExpressionData}} to generate the
##'   data and \code{interimTrial} to simulate a single survival study
##' 
##'   \code{\link{plotinterimFDR}}, \code{\link{plotinterimAPR}}, and
##'   \code{\link{plotinterimStops}} to visualize the results
##' @examples
##' ## parameters controlling the gene level expression
##' fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
##' Sigma.1 <- diag(1000)
##' Sigma.2 <- diag(1000)
##' 
##' ## parameters controlling the patient data
##' N <- 50
##' l.1.tick <- 60
##' l.2 <- 60
##' lambda <- 60
##' 
##' ## parameters controlling the study design
##' M.1 <- 2
##' M.2 <- 2
##' alpha <- 0.05
##' powerThreshold <- 0.8
##' adjustment <- "BH"
##' 
##' ## the number of studies to simulate
##' numSimRuns <- 2
##' 
##' result <- interimTrialSimulation(fc, Sigma.1, Sigma.2,
##'                                  N, l.1.tick, l.2, lambda,
##'                                  M.1, M.2, alpha, powerThreshold, adjustment,
##'                                  numSimRuns, parallel=FALSE)
##' result$resultTable
interimTrialSimulation <- function(fc, Sigma.1, Sigma.2,
                                   N, l.1.tick, l.2, lambda,
                                   M.1, M.2, alpha, powerThreshold, adjustment,
                                   numSimRuns, parallel=TRUE)
{
  if (!is.vector(fc))
    stop("the fold change vector 'fc' must be a vector")
  if(!is.numeric(fc))
    stop("the fold change vector 'fc' must be numeric")
  if(!is.matrix(Sigma.1))
    stop("the covariance matrix 'Sigma.1' must be a matrix")
  if(!is.numeric(Sigma.1))
    stop("the covariance matrix 'Sigma.1' must be numeric")
  if(!is.null(Sigma.2)) {
    if(!is.matrix(Sigma.2))
      stop("the covariance matrix 'Sigma.2' must be a matrix (or 'NULL')")
    if(!is.numeric(Sigma.2))
      stop("the covariance matrix 'Sigma.2' must be numeric (or 'NULL')")
  }

  if(!is.logical(parallel))
    stop("the parameter 'parallel' must be either 'TRUE' or 'FALSE'")

  if(length(fc) != dim(Sigma.1)[1])
    stop("dimensions of 'fc' and 'Sigma.1' do not match")
  if(!is.null(Sigma.2)) {
    if(!(length(dim(Sigma.1)) == length(dim(Sigma.2))))
      stop("dimensions of 'Sigma.1' and 'Sigma.2' do not match")
    if(!all((dim(Sigma.1) == dim(Sigma.2))))
      stop("dimensions of 'Sigma.1' and 'Sigma.2' do not match")
  }

  if(!is.numeric(N))
    stop("the sample size 'N' must be numeric")
  if(!is.vector(N))
    stop("the sample size 'N' must be a vector")
  if(!(length(N) == 1))
    stop("the sample size 'N' must be a single value")
  if(N <= 0)
    stop("the sample size 'N' must be positive")

  if(!is.numeric(l.1.tick))
    stop("the length of the recruitment period 'l.1.tick' must be numeric")
  if(!is.vector(l.1.tick))
    stop("the length of the recruitment period 'l.1.tick' must be a vector")
  if(!(length(l.1.tick) == 1))
    stop("the length of the recruitment period 'l.1.tick' must be a single value")
  if(l.1.tick <= 0)
    stop("the length of the recruitment period 'l.1.tick' must be positive")

  if(!is.numeric(lambda))
    stop("the mean survival time 'lambda' must be numeric")
  if(!is.vector(lambda))
    stop("the mean survival time 'lambda' must be a vector")
  if(!(length(lambda) == 1))
    stop("the mean survival time 'lambda' must be a single value")
  if(lambda <= 0)
    stop("the mean survival time 'lambda' must be positive")

  if(!is.numeric(l.2))
    stop("the length of the follow-up period 'l.2' must be numeric")
  if(!is.vector(l.2))
    stop("the length of the follow-up period 'l.2' must be a vector")
  if(!(length(l.2) == 1))
    stop("the length of the follow-up period 'l.2' must be a single value")
  if(l.2 <= 0)
    stop("the length of the follow-up period 'l.2' must be positive")

  if(!is.numeric(M.1))
    stop("the number of analyses during recruitment period 'M.1' must be numeric")
  if(!is.vector(M.1))
    stop("the number of analyses during recruitment period 'M.1' must be a vector")
  if(!(length(M.1) == 1))
    stop("the number of analyses during recruitment period 'M.1' must be a single value")
  if(M.1 <= 0)
    stop("the number of analyses during recruitment period 'M.1' must be positive")
  
  if(!is.numeric(M.2))
    stop("the number of analyses during follow-up period 'M.2' must be numeric")
  if(!is.vector(M.2))
    stop("the number of analyses during follow-up period 'M.2' must be a vector")
  if(!(length(M.2) == 1))
    stop("the number of analyses during follow-up period 'M.2' must be a single value")
  if(M.2 <= 0)
    stop("the number of analyses during follow-up period 'M.2' must be positive")

  if(!is.numeric(alpha))
    stop("the error level 'alpha' must be numeric")
  if(!is.vector(alpha))
    stop("the error level 'alpha' must be a vector")
  if(!(length(alpha) == 1))
    stop("the error level 'alpha' must be a single value")
  if(alpha < 0 || alpha > 1)
    stop("the error level 'alpha' must be in [0,1]")

  if(!is.numeric(powerThreshold))
    stop("the 'powerThreshold' must be numeric")
  if(!is.vector(powerThreshold))
    stop("the 'powerThreshold' must be a vector")
  if(!(length(powerThreshold) == 1))
    stop("the 'powerThreshold' must be a single value")
  if(powerThreshold < 0 || powerThreshold > 1)
    stop("the 'powerThreshold' must be in [0,1]")

  if (!(adjustment %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                          "fdr", "none")))
    stop("the adjustment method 'adjustment' must be in c(\"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\", \"none\")")

  if(!is.numeric(N))
    stop("'numSimRuns' must be numeric")
  if(!is.vector(numSimRuns))
    stop("'numSimRuns' must be a vector")
  if(!(length(numSimRuns) == 1))
    stop("'numSimRuns' must be a single value")
  if(numSimRuns <= 0)
    stop("'numSimRuns' must be positive")

  
  

  if(is.null(Sigma.2))
    Sigma.2 <- Sigma.1  

  Sigma.1 <- chol(Sigma.1)
  Sigma.2 <- chol(Sigma.2)

  if (parallel) {
    library("doMC")
    registerDoMC() # initialize parallel execution
    "%localdo%" <- get("%dopar%")
  } else {
    "%localdo%" <- get("%do%")
  }
    
  allNumbers <- foreach(r=1:numSimRuns,
                        .combine="cbind",
                        .inorder=FALSE,
                        .maxcombine=1000,
                        .packages=c("survival", "mvtnorm")) %localdo% {
    patdat <- generatePatientData(N, l.1.tick, lambda)
    gendat <- generateExpressionData(fc, Sigma.1, Sigma.2,
                       patdat$N.1, patdat$N.2, use_cholesky=TRUE)

    interimTrial(patdat,
                 gendat,
                 l.2,
                 M.1,
                 M.2,
                 alpha,
                 powerThreshold,
                 adjustment,
                 standalone=FALSE)
  }

  resultTable <- createMeanResultTable(allNumbers, M.1, M.2, numSimRuns)
  
  return(list(fc=fc,
              lambda=lambda,
              N=N,
              l.1.tick=l.1.tick,
              l.2=l.2,
              M.1=M.1,
              M.2=M.2,
              alpha=alpha,
              powerThreshold=powerThreshold,
              adjustment=adjustment,
              numSimRuns=numSimRuns,
              resultTable=resultTable))
}


##' Restructure the Simulation Output for Display
##'
##' This function restructures the output of \code{\link{interimTrial}} when run in
##' standalone mode and is also called internally by
##' \code{\link{interimTrialSimulation}}.
##'
##' @param resultTable the vector to bring in matrix form
##' @param M.1 the number of analyses during the recruitment period to be simulated
##' @param M.2 the number of analyses during the follow-up period to be simulated
##' @param numSimRuns the number of survival studies to be simulated
##' @return the parameter resultTable restructured as matrix
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
createMeanResultTable <- function(resultTable, M.1, M.2, numSimRuns=1)
{
  stoppingAnalyses <- resultTable[1,]
  stoppingAnalyses <- table(stoppingAnalyses)
  stoppingAnalysesCount <- rep(0, (M.1+M.2))
  stoppingAnalysesCount[as.integer(names(stoppingAnalyses))] <- stoppingAnalyses
  stoppingAnalysesFraction <- stoppingAnalysesCount/numSimRuns
  
  rates <- resultTable[-1,]
  if (is.vector(rates))
    meanRates <- rates
  else 
    meanRates <- rowMeans(rates)
  dim(meanRates) <- c(4,(M.1 + M.2))

  meanResultTable <- rbind(stoppingAnalysesFraction, meanRates)
  rownames(meanResultTable) <- c("Fraction of Stopped Studies",
                                 "FDR",
                                 "APR",
                                 "Estimated FDR",
                                 "Estimated APR")
  colnames(meanResultTable) <- paste("Analysis", 1:(M.1+M.2))

  return(meanResultTable)
}


##' This function visualizes the differential expression (fold changes) of the genes in
##' gendat among the groups of long survivers and short survivers.
##'
##'   The fucntion \code{generateExpressionData} has as one parameter the log fold change
##'   between the gene expression levels of the group of long surviver and
##'   the group of short survivers.  This parameter is visualized in this
##'   function with a barplot.
##'
##' @param gendat of type 'list'.  The return value of the function \code{generateExpressionData}
##' @return none
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @seealso \code{generateExpressionData}, \code{plotgeneratePatientData}
##' @examples
##' fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
##' Sigma.1 <- diag(1000)
##' Sigma.2 <- diag(1000)
##' N.1 <- 30
##' N.2 <- 30
##' 
##' gendat <- generateExpressionData(fc, Sigma.1, Sigma.2, N.1, N.2)
##' plotgendat(gendat)
plotgendat <- function(gendat)
{
  if(!is.list(gendat))
    stop("'gendat' must be the return value from the function 'generateExpressionData'")

  barplot(table(gendat$fc),
          ##cex.lab=1.5,
          ##cex.axis=1.5,
          xlab="log Fold Change",
          ylab="Frequency",
          ##cex.names=1.5,
          main="log Fold Canges Between Long and Short Survivers")
}


##' This function visualizes the patient data that is generated with
##' \code{generatePatientData}.
##'
##'   The function \code{generatePatientData} generates arrival times and survival times
##'   of patients included in the simulated survival study.  This function
##'   visualizes this data with timeline segments.  The distinction in long
##'   survivers and short survivers is color coded.
##'
##' @param patdat of type 'list'.  The return value of the function \code{generatePatientData}
##' @return none
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @seealso \code{generatePatientData}, \code{plotgenerateExpressionData}
##' @examples
##' N <- 50
##' l.1.tick <- 60
##' lambda <- 60
##' 
##' patdat <- generatePatientData(N, l.1.tick, lambda)
##' plotpatdat(patdat)
plotpatdat <- function(patdat)
{
  if(!is.list(patdat))
    stop("'patdat' must be the return value from the function 'generatePatientData'")

  N <- patdat$N.1 + patdat$N.2
  
  x0 <- patdat$arrivalTimes
  x1 <- patdat$arrivalTimes + patdat$survivalTimes
  y0 <- 1:N
  y1 <- 1:N

  arrivalOrder <- order(x0)

  x0 <- x0[arrivalOrder]
  x1 <- x1[arrivalOrder]

  survivalTimesSorted <- sort(patdat$survivalTimes)
  groupCutoff <- survivalTimesSorted[patdat$N.1]
  shortGroupIndex <- patdat$survivalTimes[arrivalOrder] <= groupCutoff
  longGroupIndex <- 1:N
  longGroupIndex <- longGroupIndex[!shortGroupIndex]

  plot(NULL,NULL,
       xlim=c(0,max(x1)),
       xlab="Time",
       ylim=c(0, N),
       ylab="Patient",
       main=paste("Patient Data"))
  
  segments(x0[shortGroupIndex],
           y0[shortGroupIndex],
           x1[shortGroupIndex],
           y1[shortGroupIndex],
           col="red",
           lwd=3)
  segments(x0[longGroupIndex],
           y0[longGroupIndex],
           x1[longGroupIndex],
           y1[longGroupIndex],
           col="green",
           lwd=3)

  legend("topright",
         legend=c("short survivers", "long survivers"),
         bg="white",
         lwd=2,
         col=c("red", "green"))
  
  abline(v=max(x0), lty="dashed")
  text(x=max(x0), y=0, labels=c("l.1"), pos=4, offset=0.5)
}


##' This function takes the result of either \code{interimTrial} or
##'  \code{interimTrialSimulation} and creates a barplot diplaying the
##'  achieved error rate (FDR) at each simulated interim analysis.
##'
##'   For each simulated interim analysis one bar is drawn that's height
##'   represents the achieved false discovery rate (FDR) at that particular
##'   interim analysis.  Additionally a dashed line is drawn at the error
##'   level alpha that was to be controlled.
##'
##' @param interimdat of type 'list'.  The return value of either \code{interimTrial} or
##'    \code{interimTrialSimulation}
##' @return none
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @seealso \code{plotinterimAPR}, \code{plotinterimStops} for other plotting
##'   functions on the same data
##' 
##'   \code{interimTrial} and \code{interimTrialSimulation} for the data generation
##' @examples
##' ## generate the data
##' fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
##' Sigma.1 <- diag(1000)
##' Sigma.2 <- diag(1000)
##' 
##' N <- 50
##' l.1.tick <- 60
##' l.2 <- 60
##' lambda <- 60
##' 
##' M.1 <- 2
##' M.2 <- 2
##' alpha <- 0.05
##' powerThreshold <- 0.8
##' adjustment <- "BH"
##' 
##' numSimRuns <- 2
##' 
##' result <- interimTrialSimulation(fc, Sigma.1, Sigma.2,
##'                                  N, l.1.tick, l.2, lambda,
##'                                  M.1, M.2, alpha, powerThreshold, adjustment,
##'                                  numSimRuns, parallel=FALSE)
##' 
##' plotinterimFDR(result)
plotinterimFDR <- function(interimdat)
{
  if(!is.list(interimdat))
    stop("'interimdat' must be the return value from the function 'interimTrial' or interimTrialSimulation")

  par(mar=c(5.1,5.1,2.1,2.1))

  plotVectorRow <- 2
  plotVector <- interimdat$resultTable[plotVectorRow,]
  plotName <- row.names(interimdat$resultTable)[plotVectorRow]

  ylimMax <- max(plotVector[!is.na(plotVector)])
  ylimMax <- max(ylimMax, 0.06)

  barplot(plotVector,
          col=rep("gray67", length(plotVector)),
          names.arg=1:(length(plotVector)),
          xlab="Interim Analysis",
          ylab=plotName,
          main=paste("Simulated", plotName),
          ylim=c(0,ylimMax))
  abline(h=interimdat$alpha, lty="dashed")
}


##' This function takes the result of either \code{interimTrial} or
##' \code{interimTrialSimulation} and creates a plot displaying the
##' achieved real and estimated average power rate (APR and eAPR) at each
##' simulated interim analysis.
##'
##'   This function plots one graph showing the achieved average power rate
##'   (APR) and one graph showing the estimated average power rate (eAPR).
##'   Both are plotted against the interim analyses such that there
##'   development over time becomes visualized.  The surivial study was
##'   stopped, when the estimated power rate (eAPR) exceeded a given
##'   threshold, which is represented by a dashed line.
##'
##' @param interimdat of type 'list'.  The return value of either \code{interimTrial} or
##'    \code{interimTrialSimulation}
##' @return none
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @seealso \code{plotinterimFDR}, \code{plotinterimStops} for other plotting
##'   functions on the same data
##' 
##'   \code{interimTrial} and \code{interimTrialSimulation} for the data generation
##' @examples
##' ## generate the data
##' fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
##' Sigma.1 <- diag(1000)
##' Sigma.2 <- diag(1000)
##' 
##' N <- 50
##' l.1.tick <- 60
##' l.2 <- 60
##' lambda <- 60
##' 
##' M.1 <- 2
##' M.2 <- 2
##' alpha <- 0.05
##' powerThreshold <- 0.8
##' adjustment <- "BH"
##' 
##' numSimRuns <- 2
##' 
##' result <- interimTrialSimulation(fc, Sigma.1, Sigma.2,
##'                                  N, l.1.tick, l.2, lambda,
##'                                  M.1, M.2, alpha, powerThreshold, adjustment,
##'                                  numSimRuns, parallel=FALSE)
##' 
##' plotinterimAPR(result)
plotinterimAPR <- function(interimdat)
{
  if(!is.list(interimdat))
    stop("'interimdat' must be the return value from the function 'interimTrial' or interimTrialSimulation")

  M <- interimdat$M.1 + interimdat$M.2

  par(mar=c(5.1,5.1,2.1,2.1))

  plotVectorRows <- c(3,5)

  plotVectors <- list()
  plotNames <- c()
  for (i in 1:(length(plotVectorRows)))
  {
    plotVectorRow <- plotVectorRows[i]
    plotVectors[[i]] <- interimdat$resultTable[plotVectorRow,]
    plotNames[i] <- row.names(interimdat$resultTable)[plotVectorRow]
  }

  ylimMax <- 1

  plot(NULL,NULL,
       xlim=c(1, M),
       xlab="Interim Analysis",
       ylim=c(0,ylimMax),
       ylab="APR",
       main=paste("Simulated", paste(plotNames, collapse=" and ")))

  legend_texts <- c()
  cols <- c()
  pch <- 1
  pchs <- c()
  lty <- 1
  ltys <- c()
  for (i in 1:(length(plotVectorRows)))
  {
    lines(1:M, plotVectors[[i]],
          type="l",
          pch=pch,
          lty=lty,
          lwd=3)
    pchs <- c(pchs, pch)
    pch <- pch + 1
    ltys <- c(ltys, lty)
    lty <- lty + 1
    legend_texts <- c(legend_texts, plotNames[i])
  }
  
  legend("bottomright",legend=legend_texts,lty=ltys, bg="white", cex=1, lwd=3)
  abline(h=interimdat$powerThreshold, lty="dashed")
}


##' This function creates a barplot showing at which interim analyses the
##' simulated survival studies were stopped.
##'
##'   For each simulated interim analysis one bar is drawn that's height
##'   represents the fraction of simulated survival studies, that were
##'   stopped at that particular interim analysis.
##'
##' @param interimdat of type 'list'.  The return value of either \code{interimTrial} or
##'    \code{interimTrialSimulation}
##' @return none
##' @export
##' @author Andreas Leha \email{andreas.leha@@med.uni-goettingen.de}
##' @seealso \code{plotinterimFDR}, \code{plotinterimAPR} for other plotting
##'   functions on the same data
##' 
##'   \code{interimTrial} and \code{interimTrialSimulation} for the data generation
##' @examples
##' ## generate the data
##' fc <- c(rep(0, 500), ceiling(rnorm(500, 0, 1)-0.5))
##' Sigma.1 <- diag(1000)
##' Sigma.2 <- diag(1000)
##' 
##' N <- 50
##' l.1.tick <- 60
##' l.2 <- 60
##' lambda <- 60
##' 
##' M.1 <- 2
##' M.2 <- 2
##' alpha <- 0.05
##' powerThreshold <- 0.8
##' adjustment <- "BH"
##' 
##' numSimRuns <- 2
##' 
##' result <- interimTrialSimulation(fc, Sigma.1, Sigma.2,
##'                                  N, l.1.tick, l.2, lambda,
##'                                  M.1, M.2, alpha, powerThreshold, adjustment,
##'                                  numSimRuns, parallel=FALSE)
##' 
##' plotinterimStops(result)
plotinterimStops <- function(interimdat)
{
  if(!is.list(interimdat))
    stop("'interimdat' must be the return value from the function 'interimTrial' or interimTrialSimulation")

  par(mar=c(5.1,5.1,2.1,2.1))

  plotVectorRow <- 1
  plotVector <- interimdat$resultTable[plotVectorRow,]
  plotName <- row.names(interimdat$resultTable)[plotVectorRow]

  ylimMax <- max(plotVector[!is.na(plotVector)])
  ylimMax <- max(ylimMax, 0.06)

  barplot(plotVector,
          col=rep("gray67", length(plotVector)),
          names.arg=1:(length(plotVector)),
          xlab="Interim Analysis",
          ylab=plotName,
          main=paste("Simulated", plotName),
          ylim=c(0,ylimMax))
}
