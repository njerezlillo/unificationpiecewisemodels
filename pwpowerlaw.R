# C's
auxiliar_pwpowerlaw <- function (p, alpha)
{
  if (any(sort(p) != p)) {
    stop(paste("p must be crecient", "\n", ""))
  }
  if (any(alpha <= 1)) {
    stop(paste("alpha must be > 1", "\n", ""))
  }
  
  aux <- c(1, rep(NA, length(p) - 1), 0)
  for (i in 2:length(p)) {
    aux[i] <- prod((p[i] / p[i - 1]) ^ (1 - alpha[i - 1]), aux[(i - 1)])
  }
  
  return (aux)
}

# density distribution function
dpwpowerlaw <- function (x, p, alpha)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar_pwpowerlaw(p, alpha)
  j <- max(which(p < x))
  d <- ((alpha[j] - 1)/p[j]) * (x/p[j])^(-alpha[j]) * C[j]

  return(d)
}

# cumulative distribution function
ppwpowerlaw <- function (x, p, alpha)
{
  C <- auxiliar_pwpowerlaw(p, alpha)
  j <- max(which(p <= x))
  u <- 1 - (x / p[j]) ^ (1 - alpha[j]) * C[j]
  
  return(u)
}

# survival function
spwpowerlaw <- function (x, p, alpha)
{
  s <- 1 - ppwpowerlaw(x, p, alpha)
  
  return(s)
}

# hazard function
hpwpowerlaw <- function (x, p, alpha)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar_pwpowerlaw(p, alpha)
  j <- max(which(p <= x))
  (alpha[j] - 1) / x
}

# quantile function
qpwpowerlaw <- function (u, p, alpha)
{
  C <- auxiliar_pwpowerlaw(p, alpha)
  aux <- u < (1 - C)
  j <- min(which(aux == TRUE)) - 1
  x <- p[j] * ((1 - u) / C[j]) ^ (1 / (1 - alpha[j]))
  
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

# sampling piece-wise power law
rpwpowerlaw <- function(n, p, alpha)
{
  q <- Vectorize(function(x) qpwpowerlaw(x, p, alpha), "x")
  while(TRUE)
  {
    x <- q(runif(n))
    nj <- n_each_interval(x, p)
    if (all(nj > 2)) break
  }
  
  return (x)
}

# Simulation under random censoring
censrand_pwpowerlaw <- function(n, p, alpha, PI)
{
  f <- function(x) x * dpwpowerlaw(x, p, alpha)
  E <- integrate(Vectorize(f), p[1], Inf)$value
  lambda <- E/PI
  q <- Vectorize(function(x) qpwpowerlaw(x, p, alpha), "x")
  
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
sum_logdata_pwpowerlaw <- function (x, p)
{
  k <- length(p)
  index <- index_each_interval(x, p)
  logdata <- vector(length = k)
  for (j in 1:k)
  {
    logdata[j] <- sum(log(x[index[[j]]]))
  }
  
  return (logdata)
}

# loglikelihood
loglik_pwpowerlaw <- function(theta, x, p)
{
  logC <- log(auxiliar_pwpowerlaw(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwpowerlaw(x, p)
  nj <- n_each_interval(x, p)
  l <-
    sum(nj * (log(theta - 1) - (1 - theta) * log(p) + logC)) - sum(theta * logdata)
  
  return (l)
}

# MLE
mle_pwpowerlaw <- function (x, p)
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
    a <- sum(log(x[index[[j]]] / p[j]))
    b <- sum(nj[(j + 1):k] * log(p[(j + 1)] / p[j]))
    alpha[j] <- 1 + nj[j] * (a + b) ^ (-1)
  }
  
  alpha[k] <- 1 + nj[k] * sum(log(x[index[[k]]] / p[k])) ^ (-1)
  
  return(alpha)
}

# profile_log_likelihood
profile_loglik_pwpowerlaw <- function(p, x)
{
  theta <- mle_pwpowerlaw(x, p)
  logC <- log(auxiliar_pwpowerlaw(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwpowerlaw(x, p)
  nj <- n_each_interval(x, p)
  l <-
    sum(nj * (log(theta - 1) - (1 - theta) * log(p) + logC)) - sum(theta * logdata)
  
  return (l)
}

# loglikelihood censored
loglik_cens_pwpowerlaw <- function(theta, x, p)
{
  logC <- log(auxiliar_pwpowerlaw(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwpowerlaw(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  l <-
    sum(dj * log(theta - 1) +
          nj * logC - nj * (1 - theta) * log(p)) +
    sum((1 - theta) * logdata) - sum(log(x$time[x$status==1]))
  
  return (l)
}

# MLE to censored data
mle_cens_pwpowerlaw <- function (x, p)
{
  k <- length(p)
  alpha <- vector(length = k)
  index <- index_each_interval(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  if (any(dj <= 2)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum(log(x$time[index[[j]]] / p[j]))
    b <- sum(nj[(j + 1):k] * log(p[(j + 1)] / p[j]))
    alpha[j] <- 1 + dj[j] * (a + b) ^ (-1)
  }
  
  alpha[k] <- 1 + dj[k] * sum(log(x$time[index[[k]]] / p[k])) ^ (-1)
  
  return(alpha)
}

# profile_log_cens_likelihood
profile_loglik_cens_pwpowerlaw <- function(p, x)
{
  theta <- mle_cens_pwpowerlaw(x, p)
  logC <- log(auxiliar_pwpowerlaw(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwpowerlaw(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  l <-
    sum(dj * log(theta - 1) +
          nj * logC - nj * (1 - theta) * log(p)) +
    sum((1 - theta) * logdata) - sum(log(x$time[x$status==1]))
  
  return (l)
}
