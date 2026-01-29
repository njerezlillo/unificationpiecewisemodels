source("pwrayleigh.R")

# Simulation Study

Table <- function(n, p, alpha) {
  R <- replicate(5000, rpwrayleigh(n, p, alpha))
  MLE <- apply(R, 2, function(x) mle_pwrayleigh(x, p))
  Nj <- apply(R, 2, function(x) n_each_interval(x, p))
  
  li <- MLE - 1.96 * alpha/sqrt(Nj)
  ls <- MLE + 1.96 * alpha/sqrt(Nj)
  
  BIAS <- round(apply(MLE - alpha, 1, mean), 3)
  VAR <- round(apply(MLE - alpha, 1, var), 3)
  CP <- round(apply(li < alpha & ls > alpha, 1, mean), 3)
  
  data.frame(BIAS, VAR, CP)
}

Table(100, p = c(0, 0.5, 1), alpha = c(1, 2, 2.3))
Table(300, p = c(0, 0.5, 1), alpha = c(1, 2, 2.3))
Table(500, p = c(0, 0.5, 1), alpha = c(1, 2, 2.3))

Table(100, p = c(0, 0.5, 1), alpha = c(2.7, 2, 2.3))
Table(300, p = c(0, 0.5, 1), alpha = c(2.7, 2, 2.3))
Table(500, p = c(0, 0.5, 1), alpha = c(2.7, 2, 2.3))

Table(100, p = c(0, 0.5, 1), alpha = c(1.5, 2, 1.5))
Table(300, p = c(0, 0.5, 1), alpha = c(1.5, 2, 1.5))
Table(500, p = c(0, 0.5, 1), alpha = c(1.5, 2, 1.5))
