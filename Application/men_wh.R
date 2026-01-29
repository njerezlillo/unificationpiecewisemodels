source("pwwh.R")
library(dplyr)
library(survival)
library(maxLik)
library(ggfortify)

x <- read.delim('./Application/skin_male.tsv', sep = '\t')
df <- data.frame(time = x$time, status = ifelse(x$censored == "true", 0, 1))
KS_fit <- survfit(Surv(time, status) ~ 1, df)
n <- nrow(df)

# Change points estimation ------------------------------------------------

target <- function(y) profile_loglik_cens_pwwh(p = c(0, y), df)

# k = 1

A_1 <- matrix(1, ncol = 1, byrow = T)
d_1 <- -1
fit_p_1 <- maxSANN(target, start = 4,
                   constraints = list(ineqA = A_1, ineqB = d_1))
fit_p_1 <- c(0, fit_p_1$estimate)
fit_alpha_1 <- mle_cens_pwwh(df, fit_p_1)

# k = 2

A_2 <- matrix(c(1, 0, -1, 1, 0, -1), ncol = 2, byrow = T)
d_2 <- c(-0.5, -1, 20)
fit_p_2 <- maxSANN(target, start = c(3, 6),
                   constraints = list(ineqA = A_2, ineqB = d_2))
fit_p_2 <- c(0, fit_p_2$estimate)
fit_alpha_2 <- mle_cens_pwwh(df, fit_p_2)

# k = 3

A_3 <- matrix(c(1, 0, 0, -1, 1, 0, 0, -1, 1), ncol = 3, byrow = T)
d_3 <- c(-0.5, -1, -1)
fit_p_3 <- maxSANN(target, start = c(1, 4, 6),
                   constraints = list(ineqA = A_3, ineqB = d_3))
fit_p_3 <- c(0, fit_p_3$estimate)
fit_alpha_3 <- mle_cens_pwwh(df, fit_p_3)

# AIC ---------------------------------------------------------------------

AIC_1 <- 2 * 3 - 2 * loglik_cens_pwwh(fit_alpha_1, df, fit_p_1)
AIC_2 <- 2 * 5 - 2 * loglik_cens_pwwh(fit_alpha_2, df, fit_p_2)
AIC_3 <- 2 * 7 - 2 * loglik_cens_pwwh(fit_alpha_3, df, fit_p_3)

# BIC ---------------------------------------------------------------------

BIC_1 <- log(n) * 3 - 2 * loglik_cens_pwwh(fit_alpha_1, df, fit_p_1)
BIC_2 <- log(n) * 5 - 2 * loglik_cens_pwwh(fit_alpha_2, df, fit_p_2)
BIC_3 <- log(n) * 7 - 2 * loglik_cens_pwwh(fit_alpha_3, df, fit_p_3)

# Output ------------------------------------------------------------------

gof <- data.frame(
  "k" = 1:3,
  AIC = c(AIC_1, AIC_2, AIC_3),
  BIC = c(BIC_1, BIC_2, BIC_3)
)

gof

# print(xtable(gof), include.rownames = F)

# Sequential testing ------------------------------------------------------

significance <- 0.01
var <- (fit_alpha_1)^2/n_each_interval(df$time[df$status == 1], fit_p_1)
w <- (fit_alpha_1[1] - fit_alpha_1[2])^2/sum(var[1:2])
pchisq(w, 1, lower.tail = F) < significance

var <- (fit_alpha_2)^2/n_each_interval(df$time[df$status == 1], fit_p_2)
w1 <- (fit_alpha_2[1] - fit_alpha_2[2])^2/sum(var[1:2])
w2 <- (fit_alpha_2[2] - fit_alpha_2[3])^2/sum(var[2:3])
pchisq(w1, 1, lower.tail = F) < significance/2
pchisq(w2, 1, lower.tail = F) < significance/2

var <- (fit_alpha_3)^2/n_each_interval(df$time[df$status == 1], fit_p_3)
w1 <- (fit_alpha_3[1] - fit_alpha_3[2])^2/sum(var[1:2])
w2 <- (fit_alpha_3[2] - fit_alpha_3[3])^2/sum(var[2:3])
w3 <- (fit_alpha_3[3] - fit_alpha_3[4])^2/sum(var[3:4])
pchisq(w1, 1, lower.tail = F) < significance/4
pchisq(w2, 1, lower.tail = F) < significance/4
pchisq(w3, 1, lower.tail = F) < significance/4
