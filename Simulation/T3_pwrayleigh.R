source("pwrayleigh.R")

# Part 1 ------------------------------------------------------------------

p = c(0, 0.5, 1)
PI = .3
alpha = c(1, 2, 2.3)

R1_100 <- replicate(5000, censrand_pwrayleigh(100, p, alpha, PI), simplify = FALSE)
R1_300 <- replicate(5000, censrand_pwrayleigh(300, p, alpha, PI), simplify = FALSE)
R1_500 <- replicate(5000, censrand_pwrayleigh(500, p, alpha, PI), simplify = FALSE)

alpha = c(2.7, 2, 2.3)
R2_100 <- replicate(5000, censrand_pwrayleigh(100, p, alpha, PI), simplify = FALSE)
R2_300 <- replicate(5000, censrand_pwrayleigh(300, p, alpha, PI), simplify = FALSE)
R2_500 <- replicate(5000, censrand_pwrayleigh(500, p, alpha, PI), simplify = FALSE)

alpha = c(1.5, 2, 1.5)
R3_100 <- replicate(5000, censrand_pwrayleigh(100, p, alpha, PI), simplify = FALSE)
R3_300 <- replicate(5000, censrand_pwrayleigh(300, p, alpha, PI), simplify = FALSE)
R3_500 <- replicate(5000, censrand_pwrayleigh(500, p, alpha, PI), simplify = FALSE)

#save(R1_100,R1_300,R1_500, file="./Simulation/censrand_st1_rayleigh.RData")
#save(R2_100,R2_300,R2_500, file="./Simulation/censrand_st2_rayleigh.RData")
#save(R3_100,R3_300,R3_500, file="./Simulation/censrand_st3_rayleigh.RData")

# Part 2 ------------------------------------------------------------------

Table <- function(n, p, alpha,R) {
  MLE <- do.call(rbind, lapply(R, function(x) mle_cens_pwrayleigh(x, p)))
  Nj <- do.call(rbind, lapply(R, function(x) n_each_interval(x$time[x$status==1], p)))
  
  li <- MLE - 1.96 * MLE / sqrt(Nj)
  ls <- MLE + 1.96 * MLE / sqrt(Nj)
  
  BIAS <- round(apply(t(MLE) - alpha, 1, mean), 3)
  VAR <- round(apply(t(MLE) - alpha, 1, var), 3)
  CP <- round(apply(t(li) < alpha & t(ls) > alpha, 1, mean), 3)
  
  data.frame(BIAS, VAR, CP)
}

load("./Simulation/censrand_st1_rayleigh.RData")
load("./Simulation/censrand_st2_rayleigh.RData")
load("./Simulation/censrand_st3_rayleigh.RData")

# Escenario 1
Table(100, p = c(0, 0.5, 1), alpha = c(1, 2, 2.3), R1_100) 
Table(300, p = c(0, 0.5, 1), alpha = c(1, 2, 2.3), R1_300) 
Table(500, p = c(0, 0.5, 1), alpha = c(1, 2, 2.3), R1_500)

# Escenario 2
Table(100, p = c(0, 0.5, 1), alpha = c(2.7, 2, 2.3), R2_100) 
Table(300, p = c(0, 0.5, 1), alpha = c(2.7, 2, 2.3), R2_300) 
Table(500, p = c(0, 0.5, 1), alpha = c(2.7, 2, 2.3), R2_500)

# Escenario 3
Table(100, p = c(0, 0.5, 1), alpha = c(1.5, 2, 1.5), R3_100) 
Table(300, p = c(0, 0.5, 1), alpha = c(1.5, 2, 1.5), R3_300) 
Table(500, p = c(0, 0.5, 1), alpha = c(1.5, 2, 1.5), R3_500)
