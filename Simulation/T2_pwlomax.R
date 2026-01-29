source("pwlomax.R")

# Simulation Study

Table <- function(n, p, alpha) {
  R <- replicate(5000, rpwlomax(n, p, alpha))
  MLE <- apply(R, 2, function(x) mle_pwlomax(x, p))
  Nj <- apply(R, 2, function(x) n_each_interval(x, p))
  
  li <- MLE - 1.96 * alpha/sqrt(Nj)
  ls <- MLE + 1.96 * alpha/sqrt(Nj)
  
  BIAS <- round(apply(MLE - alpha, 1, mean), 3)
  RMSE <- round(apply(MLE - alpha, 1, function(x) sqrt(mean(x^2))), 3)
  CP <- round(apply(li < alpha & ls > alpha, 1, mean), 3)
  
  data.frame(BIAS, RMSE, CP)
}

Table(100, p = c(0, 0.7, 1.4), alpha = c(0.8, 0.6, 1.3))
Table(300, p = c(0, 0.7, 1.4), alpha = c(0.8, 0.6, 1.3))
Table(500, p = c(0, 0.7, 1.4), alpha = c(0.8, 0.6, 1.3))

Table(100, p = c(0, 0.7, 1.4), alpha = c(1.2, 0.7, 1.1))
Table(300, p = c(0, 0.7, 1.4), alpha = c(1.2, 0.7, 1.1))
Table(500, p = c(0, 0.7, 1.4), alpha = c(1.2, 0.7, 1.1))

Table(100, p = c(0, 0.7, 1.4), alpha = c(0.8, 1.5, 1.2))
Table(300, p = c(0, 0.7, 1.4), alpha = c(0.8, 1.5, 1.2))
Table(500, p = c(0, 0.7, 1.4), alpha = c(0.8, 1.5, 1.2))
