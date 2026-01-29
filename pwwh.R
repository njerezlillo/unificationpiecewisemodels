# C's
auxiliar_pwwh <- function (p, beta)
{
  if (any(sort(p) != p)) {
    stop(paste("p must be crecient", "\n", ""))
  }
  if (any(beta < 0)) {
    stop(paste("beta must be > 0", "\n", ""))
  }
  
  aux <- c(1, rep(NA, length(p) - 1), 0)
  for (i in 2:length(p)) {
    aux[i] <- prod(exp((-1/beta)[i - 1]*(p[i] - p[i - 1])^3), aux[(i - 1)])
  }
  
  return (aux)
}

# density distribution function
dpwwh <- function (x, p, beta)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar_pwwh(p, beta)
  j <- max(which(p < x))
  d <- 3*(x - p[j])^2*(1/beta[j])* exp(-1/beta[j] * (x - p[j])^3) * C[j]
  
  return(d)
}

# cumulative distribution function
ppwwh <- function (x, p, beta)
{
  C <- auxiliar_pwwh(p, beta)
  j <- max(which(p <= x))
  u <- 1 - exp(-(1/beta[j]) * (x - p[j])^3) * C[j]
  
  return(u)
}

# survival function
spwwh <- function (x, p, beta)
{
  C <- auxiliar_pwwh(p, beta)
  j <- max(which(p <= x))
  s <- exp(-1/beta[j] * (x - p[j])^3) * C[j]
  
  return(s)
}

# hazard function
hpwwh <- function (x, p, beta)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar_pwwh(p, beta)
  j <- max(which(p <= x))
  h<- (3/beta[j])*(x - p[j])^2
  
}

# quantile function
qpwwh <- function (u, p, beta)
{
  C <- auxiliar_pwwh(p, beta)
  aux <- u < (1 - C)
  j <- min(which(aux == TRUE)) - 1
  x<- p[j] + (-beta[j]*log((1 - u) / C[j])) ^ (1/3)
  return(x)
}


# identifying the observations in each interval
index_each_interval <- function (x, p)
{
  k <- length(p)
  nj <- vector("list", length = k)
  for (j in 1:k)
  {
    if (j < k)
    {
      nj[[j]] <- which(x >= p[j] & x < p[j + 1])
    } else {
      nj[[j]] <- which(x >= p[j])
    }
  }
  
  return (nj)
}

# counting the observations in each interval
n_each_interval <- function (x, p)
{
  nj <- lengths(index_each_interval(x, p))
  return (nj)
}


# sampling piecewise wilson-hilferty 
rpwwh <- function(n, p, beta)
{
  q <- Vectorize(function(x) qpwwh(x, p, beta), "x")
  while(TRUE)
  {
    x <- q(runif(n))
    nj <- n_each_interval(x, p)
    if (all(nj > 2)) break
  }
  
  return (x)
}

# Simulation under random censoring
censrand_pwwh <- function(n, p, beta, PI)
{
  f <- function(x) x * dpwwh(x, p, beta)
  E <- integrate(Vectorize(f), p[1], Inf)$value
  lambda <- E/PI
  delta <- vector(length = n) ####
  C <- auxiliar_pwwh(p, beta) ####
  q <- Vectorize(function(x) qpwwh(x, p, beta), "x")
  
  while(TRUE)
  {
    x <- q(runif(n))
    v <- runif(n, p[1], lambda)
    delta <- (x < v)
    x[!delta] <- v[!delta]
    
    out <- data.frame(time = x, status = as.numeric(delta))
    dj <- n_each_interval(out$time[out$status == 1], p)
    if (all(dj > 2)) break
  }
  
  return (out)
}

# sum(logdata) in each I[j]
sum_logdata_pwwh <- function (x, p)
{
  k <- length(p)
  index <- index_each_interval(x, p)
  logdata1 <- vector(length = k)
  logdata2<- vector(length = k)
  for (j in 1:k)
  {
    logdata1[j] <- sum((x[index[[j]]] - p[j])^3)
    logdata2[j] <- sum(log((x[index[[j]]] - p[j])^2))
  }
  logdata=rbind(logdata1,logdata2)
  return (logdata)
}

# log_likelihood
loglik_pwwh <- function(theta, x, p)
{
  logC <- log(auxiliar_pwwh(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwwh(x, p)
  nj <- n_each_interval(x, p)
  l <- sum(nj * (log(3/theta) + logC)) + sum(logdata[2,])- sum(logdata[1,]/(theta))
  return (l)
}

# MLE
mle_pwwh <- function (x, p)
{
  k <- length(p)
  beta <- vector(length = k)
  index <- index_each_interval(x, p)
  nj <- n_each_interval(x, p)
  if (any(nj == 0)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum((x[index[[j]]] - p[j])^3)
    b <- sum(nj[(j + 1):k] * (p[(j + 1)] - p[j])^3)
    beta[j] <- (a + b)/(nj[j])
  }
  
  beta[k] <- sum((x[index[[k]]] - p[k])^3)/(nj[k])
  
  return(beta)
}

# profile_log_likelihood
profile_loglik_pwwh <- function(p, x)
{
  theta <- mle_pwwh(x, p)
  logC <- log(auxiliar_pwwh(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwwh(x, p)
  nj <- n_each_interval(x, p)
  l <- sum(nj * (log(3/theta) + logC)) + sum(logdata[2,])- sum(logdata[1,]/(theta))
  return (l)
}

# loglikelihood censored
loglik_cens_pwwh <- function(theta, x, p)
{
  logC <- log(auxiliar_pwwh(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwwh(x$time, p)
  logdata[2,] <- sum_logdata_pwwh(x$time[x$status == 1], p)[2,]
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  l <-sum(dj * log(3/theta) + nj * logC) + sum(logdata[2,])- sum(logdata[1,]/(theta))
  return (l)
}

# MLE to censored data
mle_cens_pwwh <- function (x, p)
{
  k <- length(p)
  beta <- vector(length = k)
  index <- index_each_interval(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  if (any(dj == 0)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum((x$time[index[[j]]] - p[j])^3)
    b <- sum(nj[(j + 1):k] * (p[(j + 1)] - p[j])^3)
    beta[j] <-(a + b)/dj[j]  
  }
  
  beta[k] <- sum((x$time[index[[k]]] - p[k])^3)/dj[k]
  
  return(beta)
}

# profile_log_cens_likelihood
profile_loglik_cens_pwwh <- function(p, x)
{
  theta <- mle_cens_pwwh(x, p)
  logC <- log(auxiliar_pwwh(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwwh(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  l <-sum(dj * log(3/theta) + nj * logC) + sum(logdata[2,])- sum(logdata[1,]/(theta))
  return (l)
}
