source("pwrayleigh.R")
source("pwexp.R")
library(dplyr)
library(survival)
library(maxLik)
library(ggfortify)
library(flexsurv)
library(survminer)

x <- read.delim('./Application/skin_female.tsv', sep = '\t')
df <- data.frame(time = x$time, status = ifelse(x$censored == "true", 0, 1))
KS_fit <- survfit(Surv(time, status) ~ 1, df)
n <- length(df$time)

# Change points estimation ------------------------------------------------

target <- function(y) profile_loglik_cens_pwrayleigh(p = c(0, y), df)

# k = 1

A_1 <- matrix(1, ncol = 1, byrow = T)
d_1 <- -1
fit_p_1 <- maxSANN(target, start = 4,
                   constraints = list(ineqA = A_1, ineqB = d_1))
fit_p_1 <- c(0, fit_p_1$estimate)
fit_alpha_1 <- mle_cens_pwrayleigh(df, fit_p_1)

# k = 2

A_2 <- matrix(c(1, 0, -1, 1, 0, -1), ncol = 2, byrow = T)
d_2 <- c(-0.5, -1, 20)
fit_p_2 <- maxSANN(target, start = c(3, 6),
                   constraints = list(ineqA = A_2, ineqB = d_2))
fit_p_2 <- c(0, fit_p_2$estimate)
fit_alpha_2 <- mle_cens_pwrayleigh(df, fit_p_2)

# k = 3

A_3 <- matrix(c(1, 0, 0, -1, 1, 0, 0, -1, 1), ncol = 3, byrow = T)
d_3 <- c(-0.5, -1, -1)
fit_p_3 <- maxSANN(target, start = c(1, 4, 6),
                   constraints = list(ineqA = A_3, ineqB = d_3))
fit_p_3 <- c(0, fit_p_3$estimate)
fit_alpha_3 <- mle_cens_pwrayleigh(df, fit_p_3)

# AIC ---------------------------------------------------------------------

AIC_1 <- -2 * loglik_cens_pwrayleigh(fit_alpha_1, df, fit_p_1) + 2 * 3
AIC_2 <- -2 * loglik_cens_pwrayleigh(fit_alpha_2, df, fit_p_2) + 2 * 5
AIC_3 <- -2 * loglik_cens_pwrayleigh(fit_alpha_3, df, fit_p_3) + 2 * 7

# BIC ---------------------------------------------------------------------

BIC_1 <- -2 * loglik_cens_pwrayleigh(fit_alpha_1, df, fit_p_1) + log(n) * 3
BIC_2 <- -2 * loglik_cens_pwrayleigh(fit_alpha_2, df, fit_p_2) +  log(n) * 5
BIC_3 <- -2 * loglik_cens_pwrayleigh(fit_alpha_3, df, fit_p_3) + log(n) * 7

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

# Ci ----------------------------------------------------------------------

se <- fit_alpha_2/sqrt(n_each_interval(df$time[df$status == 1], fit_p_2))
ic_l <- round(fit_alpha_2 - 1.96 * se, 2)
ic_u <- round(fit_alpha_2 + 1.96 * se, 2)

fit <- data.frame(j = 1:3,
                  "estimation" = fit_alpha_2,
                  "CI" = paste0("(", ic_l, ", ", ic_u, ")"))
fit

# Fig ---------------------------------------------------------------------

S_ray_x <- function(x) { 
  spwrayleigh(x, c(0, 2.396783, 5.566070), c(9.36337, 11.53082, 178.31230))
}
S_ray_x <- Vectorize(S_ray_x)

p_df_x <- data.frame(p = c(0, 2.396783, 5.566070),
                     e = S_ray_x(c(0, 2.396783, 5.566070)))

S_exp_x <- function(x) { 
  spwexp(x, c(0, 0.7214562, 5.5583778), c(0.01664315, 0.14554270, 0.03739037))
}
S_exp_x <- Vectorize(S_exp_x)

g_km_x <- 
  autoplot(KS_fit, conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.1) + theme_bw() +
  stat_function(fun=S_ray_x, col = "dodgerblue3", size = 0.7) + 
  #stat_function(fun=S_exp_x, col = "red", size = 0.7, lty = 2) + 
  geom_point(aes(x = p, y = e), col = "dodgerblue3", data = p_df_x, size = 1.2) +
  labs(x = "Time", y = "Survival") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15))
g_km_x

ggsave(plot = g_km_x, "FIT.pdf", height = 6, width = 10)

ei_x <- -log(S_ray_x(df$time)) 
km_ei_x <- survfit(Surv(ei_x, df$status) ~ 1)

g_cs_x <-
  ggplot() + aes(x = km_ei_x$surv, y = exp(-km_ei_x$time)) + 
  geom_point() + theme_bw() + 
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1.2) +
  labs(x = "S(ei): Kaplan-Meier", y = "S(ei):standard exponential") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15))

ggsave(plot = g_cs_x, "CS.pdf", height = 6, width = 10)
