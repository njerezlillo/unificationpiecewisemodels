source("pwwh.R")

# Simulation Study

Table <- function(n, p, alpha) {
  R <- replicate(5000, rpwwh(n, p, alpha))
  MLE <- apply(R, 2, function(x) mle_pwwh(x, p))
  Nj <- apply(R, 2, function(x) n_each_interval(x, p))
  
  li <- MLE - 1.96 * alpha/sqrt(Nj)
  ls <- MLE + 1.96 * alpha/sqrt(Nj)
  
  BIAS <- round(apply(MLE - alpha, 1, mean), 3)
  VAR <- round(apply(MLE - alpha, 1, var), 3)
  CP <- round(apply(li < alpha & ls > alpha, 1, mean), 3)
  
  data.frame(BIAS, VAR, CP)
}

Table(100, p = c(0, 1, 2), alpha = c(4, 7, 8))
Table(300, p = c(0, 1, 2), alpha = c(4, 7, 8))
Table(500, p = c(0, 1, 2), alpha = c(4, 7, 8))

Table(100, p = c(0, 1, 2), alpha = c(10, 8, 6))
Table(300, p = c(0, 1, 2), alpha = c(10, 8, 6))
Table(500, p = c(0, 1, 2), alpha = c(10, 8, 6))

Table(100, p = c(0, 1, 2), alpha = c(6, 8, 7))
Table(300, p = c(0, 1, 2), alpha = c(6, 8, 7))
Table(500, p = c(0, 1, 2), alpha = c(6, 8, 7))

