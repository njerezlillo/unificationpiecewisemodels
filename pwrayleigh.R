# C's
auxiliar_pwrayleigh<- function (p, alpha)
{
  if (any(sort(p) != p)) {
    stop(paste("p must be crecient", "\n", ""))
  }
  if (any(alpha <= 0)) {
    stop(paste("alpha must be > 0", "\n", ""))
  }
  
  aux <- c(1, rep(NA, length(p) - 1), 0)
  for (i in 2:length(p)) {
    aux[i] <- prod(exp(-(p[i]-p[i - 1]) ^ 2/(2*alpha[i - 1])), aux[(i - 1)])
  }
  
  return (aux)
}

# density distribution function
dpwrayleigh <- function (x, p, alpha)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar_pwrayleigh(p, alpha)
  j <- max(which(p < x))
  d <- ( (x-p[j])/alpha[j])*exp(-(x-p[j])^2 / (2*alpha[j])) * C[j]
  
  return(d)
}

# cumulative distribution function
ppwrayleigh <- function (x, p, alpha)
{
  C <- auxiliar_pwrayleigh(p, alpha)
  j <- max(which(p <= x))
  u <- 1 - exp(-(x-p[j])^2 / (2*alpha[j])) * C[j]
  
  return(u)
}

# survival function
spwrayleigh <- function (x, p, alpha)
{
  s <- 1 - ppwrayleigh(x, p, alpha)
  
  return(s)
}

# hazard function
hpwrayleigh <- function (x, p, alpha)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar_pwrayleigh(p, alpha)
  j <- max(which(p <= x))
  (x-p[j])/alpha[j]
}

# quantile function
qpwrayleigh <- function (u, p, alpha)
{
  C <- auxiliar_pwrayleigh(p, alpha)
  aux <- u < (1 - C)
  j <- min(which(aux == TRUE)) - 1
  x <- p[j] + (-2*alpha[j]*log((1 - u) / C[j])) ^ (1/2)
  
  return (x)
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

# sampling piece-wise rayleigh
rpwrayleigh <- function(n, p, alpha)
{
  q <- Vectorize(function(x) qpwrayleigh(x, p, alpha), "x")
  while(TRUE)
  {
    x <- q(runif(n))
    nj <- n_each_interval(x, p)
    if (all(nj > 2)) break
  }
  
  return (x)
}

# Simulation under random censoring
censrand_pwrayleigh <- function(n, p, alpha, PI)
{
  f <- function(x) x * dpwrayleigh(x, p, alpha)
  E <- integrate(Vectorize(f), p[1], Inf)$value
  lambda <- E/PI
  delta <- vector(length = n)
  C <- auxiliar_pwrayleigh(p, alpha)
  q <- Vectorize(function(x) qpwrayleigh(x, p, alpha), "x")
  
  while(TRUE)
  {
    x <- q(runif(n))
    v <- runif(n, 0, lambda)
    delta <- (x < v)
    x[!delta] <- v[!delta]
    
    out <- data.frame(time = x, status = as.numeric(delta))
    dj <- n_each_interval(out$time[out$status == 1], p)
    if (all(dj > 2)) break
  }
  
  return (out)
}

# sum(logdata) in each I[j]
sum_logdata_pwrayleigh<- function (x, p)
{
  k <- length(p)
  index <- index_each_interval(x, p)
  logdata1 <- vector(length = k)
  logdata2<- vector(length = k)
  for (j in 1:k)
  {
    logdata1[j] <- sum((x[index[[j]]] - p[j])^2)
    logdata2[j] <- sum(log(x[index[[j]]] - p[j]))
  }
  logdata=rbind(logdata1,logdata2)
  return (logdata)
}

# log_likelihood
loglik_pwrayleigh <- function(theta, x, p)
{
  logC <- log(auxiliar_pwrayleigh(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwrayleigh(x, p)
  nj <- n_each_interval(x, p)
  l <-
    sum(nj * (-log(theta) + logC)) + sum(logdata[2,])- sum(logdata[1,]/(2*theta))
  
  return (l)
}

# MLE
mle_pwrayleigh <- function (x, p)
{
  k <- length(p)
  alpha <- vector(length = k)
  index <- index_each_interval(x, p)
  nj <- n_each_interval(x, p)
  if (any(nj == 0)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum((x[index[[j]]] - p[j])^2)
    b <- sum(nj[(j + 1):k] * (p[(j + 1)] - p[j])^2)
    alpha[j] <- (a + b)/(2*nj[j])
  }
  
  alpha[k] <- sum((x[index[[k]]] - p[k])^2)/(2*nj[k])
  
  return(alpha)
}

# profile_log_likelihood
profile_loglik_pwrayleigh <- function(p, x)
{
  theta <- mle_pwrayleigh(x, p)
  logC <- log(auxiliar_pwrayleigh(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwrayleigh(x, p)
  nj <- n_each_interval(x, p)
  l <-
    sum(nj * (-log(theta) + logC)) + sum(logdata[2,])- sum(logdata[1,]/(2*theta))
  
  return (l)
}

## loglikelihood censored
loglik_cens_pwrayleigh <- function(theta, x, p)
{
  logC <- log(auxiliar_pwrayleigh(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwrayleigh(x$time, p)
  logdata[2,] <- sum_logdata_pwrayleigh(x$time[x$status == 1], p)[2,]
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  l <-
    sum(dj * (-log(theta) ) + nj*(logC)) + sum(logdata[2,])- sum(logdata[1,]/(2*theta))
  
  return (l)
}

# MLE to censored data
mle_cens_pwrayleigh <- function (x, p)
{
  k <- length(p)
  alpha <- vector(length = k)
  index <- index_each_interval(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  if (any(dj == 0)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum((x$time[index[[j]]] - p[j])^2)
    b <- sum(nj[(j + 1):k] * (p[(j + 1)] - p[j])^2)
    alpha[j] <- (a + b)/(2*dj[j])
  }
  
  alpha[k] <- sum((x$time[index[[k]]] - p[k])^2)/(2*dj[k])
  
  return(alpha)
}

# profile_log_cens_likelihood
profile_loglik_cens_pwrayleigh <- function(p, x)
{
  theta <- mle_cens_pwrayleigh(x, p)
  logC <- log(auxiliar_pwrayleigh(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwrayleigh(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  l <-
    sum(dj * (-log(theta) ) + nj*(logC)) + sum(logdata[2,])- sum(logdata[1,]/(2*theta))
  
  return (l)
}
