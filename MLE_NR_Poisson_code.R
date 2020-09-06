# MAXIMUM LIKELIHOOD ESTIMATION: NEWTON-RAPHSON ALGORITHM
# POISSON DISTRIBUTION

MLE_NR_Poisson <- function(y, maxiter = 100, epsilon = 0.000000001, stop_criteria = 10000){

  n <- length(y)
  sumy <- sum(y)

  U <- function(n, sumy, theta){
    (1/theta)*sumy - n
  }

  Uprime <- function(n, sumy, theta){
    -(1/(theta^2))*sumy
  }

  # first guess:
  theta <- min(y)
  Estimator <- matrix(NA, ncol = 1)
  Estimator[1, 1] <- theta

  iter <- 1
  while( (stop_criteria > epsilon) & (iter <= maxiter) ){

    num <- U(n, sumy, Estimator[iter, 1])
    den <- Uprime(n, sumy, Estimator[iter, 1])

    UPD <- -(num/den)
    Estimator_iter <- as.matrix(Estimator[iter, 1]) + UPD
    Estimator <- rbind(Estimator, Estimator_iter)
    stop_criteria <- UPD^2
    iter <- iter + 1
  }

  Likelihood <- matrix(NA, ncol = 1)
  LogLikelihood <- matrix(NA, ncol = 1)

  for(i in 1:length(Estimator)){
    Likelihood[i] <- prod( ( (Estimator[i,1]^y)*exp(-Estimator[i,1]) )/(factorial(y)) )
  }

  for(i in 1:length(Estimator)){
    LogLikelihood[i] <- sum( y*log(Estimator[i,1]) - Estimator[i,1] - log(factorial(y)) )
  }

  results <- cbind(format(Estimator, nsmall = 6), Likelihood, LogLikelihood)
  colnames(results) <- c("ML Estimator", "Likelihood", "Log-Likelihood")
  return(results)

}

