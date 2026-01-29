source("pwpowerlaw.R")

# Simulation Study

Table <- function(n, p, alpha) {
  R <- replicate(5000, rpwpowerlaw(n, p, alpha))
  MLE <- apply(R, 2, function(x) mle_pwpowerlaw(x, p))
  Nj <- apply(R, 2, function(x) n_each_interval(x, p))
  
  li <- MLE - 1.96 * (alpha - 1)/sqrt(Nj)
  ls <- MLE + 1.96 * (alpha - 1)/sqrt(Nj)
  
  BIAS <- round(apply(MLE - alpha, 1, mean), 3)
  VAR <- round(apply(MLE - alpha, 1, var), 3)
  CP <- round(apply(li < alpha & ls > alpha, 1, mean), 3)
  
  data.frame(BIAS, VAR, CP)
}

Table(100, p = c(1, 1.5, 2), alpha = c(1.8, 2, 3))
Table(300, p = c(1, 1.5, 2), alpha = c(1.8, 2, 3))
Table(500, p = c(1, 1.5, 2), alpha = c(1.8, 2, 3))

Table(100, p = c(1, 1.5, 2), alpha = c(3, 2, 2.2))
Table(300, p = c(1, 1.5, 2), alpha = c(3, 2, 2.2))
Table(500, p = c(1, 1.5, 2), alpha = c(3, 2, 2.2))

Table(100, p = c(1, 1.5, 2), alpha = c(2, 3, 2.5))
Table(300, p = c(1, 1.5, 2), alpha = c(2, 3, 2.5))
Table(500, p = c(1, 1.5, 2), alpha = c(2, 3, 2.5))
