source("pwexp.R")
source("pwrayleigh.R")
source("pwpowerlaw.R")
library(dplyr)
library(survival)
library(maxLik)
library(ggfortify)

x <- read.delim('./Application/skin_female.tsv', sep = '\t')
y <- read.delim('./Application/skin_male.tsv', sep = '\t')
df_x <- data.frame(time = x$time, status = ifelse(x$censored == "true", 0, 1))
df_y <- data.frame(time = y$time, status = ifelse(y$censored == "true", 0, 1))

KS_fit_x <- survfit(Surv(time, status) ~ 1, df_x)
KS_fit_y <- survfit(Surv(time, status) ~ 1, df_y)

# Fig 1 -------------------------------------------------------------------

S_ray_x <- function(x) { 
  spwrayleigh(x, c(0, 2.396783, 5.566070), c(9.36337, 11.53082, 178.31230))
}
S_ray_x <- Vectorize(S_ray_x)

S_exp_x <- function(x) {
  spwexp(x, c(0, 0.7214562, 5.5583778), c(0.01664315, 0.14554270, 0.03739037))
}
S_exp_x <- Vectorize(S_exp_x)

p_df_x <- data.frame(p = c(0, 2.396783, 5.566070),
                   e = S_ray_x(c(0, 2.396783, 5.566070)))

g_km_x <- 
  autoplot(KS_fit_x, conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.1) + theme_bw() +
  stat_function(fun=S_ray_x, col = "dodgerblue3", size = 0.7) + 
  #stat_function(fun=S_exp_x, col = "red", size = 0.7) + 
  geom_point(aes(x = p, y = e), col = "dodgerblue3", data = p_df_x, size = 1.2) +
  labs(x = "time", y = "survival") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15))

S_exp_y <- function(x) {
  spwexp(x, c(0, 0.719483, 2.484443), c(0.04026234, 0.14855872, 0.09460333))
}
S_exp_y <- Vectorize(S_exp_y)

p_df_y <- data.frame(p = c(0, 0.719483, 2.484443),
                     e = S_exp_y(c(0, 0.719483, 2.484443)))

g_km_y <- autoplot(KS_fit_y, conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.1) + theme_bw() +
  stat_function(fun = S_exp_y, col = "dodgerblue3", size = 0.7) +
  geom_point(aes(x = p, y = e), col = "dodgerblue3", data = p_df_y, size = 1.2) +
  labs(x = "time", y = "survival") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15))

g_km <- gridExtra::grid.arrange(g_km_x, g_km_y, ncol = 2)

#ggsave(plot = g_km, "./Application/Fig1.pdf", height = 6, width = 15)
ggsave(plot = g_km, "./Application/Fig1.eps", height = 6, width = 15)

# Fig 2 -------------------------------------------------------------------

ei_x <- -log(S_ray_x(df_x$time)) 
km_ei_x <- survfit(Surv(ei_x, df_x$status) ~ 1)

ei_y <- -log(S_exp_y(df_y$time)) 
km_ei_y <- survfit(Surv(ei_y, df_y$status) ~ 1)

g_cs_x <-
  ggplot() + aes(x = km_ei_x$surv, y = exp(-km_ei_x$time)) + 
  geom_point() + theme_bw() + 
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1.2) +
  labs(x = "S(ei): Kaplan-Meier", y = "S(ei):standard exponential") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15))

g_cs_y <-
  ggplot() + aes(x = km_ei_y$surv, y = exp(-km_ei_y$time)) + 
  geom_point() + theme_bw() + 
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1.2) +
  labs(x = "S(ei): Kaplan-Meier", y = "S(ei):standard exponential") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15))

g_cs <- gridExtra::grid.arrange(g_cs_x, g_cs_y, ncol = 2)

#ggsave(plot = g_cs, "./Application/Fig2.pdf", height = 6, width = 15)
ggsave(plot = g_cs, "./Application/Fig2.eps", height = 6, width = 15)

# Fig 3 -------------------------------------------------------------------

mef_ray_x <- function (u) {
  aux <- Vectorize(function(t) spwrayleigh(t, c(0, 2.396783, 5.566070), c(9.36337, 11.53082, 178.31230)))
  integrate(aux, u, Inf)$value / spwrayleigh(u, c(0, 2.396783, 5.566070), c(9.36337, 11.53082, 178.31230))
}

mef_ray_x <- Vectorize(mef_ray_x)

mef_exp_y <- function (u) {
  aux <- Vectorize(function(t) spwexp(t, c(0, 0.719483, 2.484443), c(0.04026234, 0.14855872, 0.09460333)))
  integrate(aux, u, Inf)$value / spwexp(u, c(0, 0.719483, 2.484443), c(0.04026234, 0.14855872, 0.09460333))
}

mef_exp_y <- Vectorize(mef_exp_y)

g_mef_x <- ggplot() + 
  stat_function(fun = mef_ray_x, col = "dodgerblue3", size = 0.7, xlim = c(0, 15)) +
  theme_bw() + scale_y_continuous(limit = c(9, 17), breaks = seq(9, 17, 2)) +
  labs(x = "time", y = "mean residual life") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15))

g_mef_y <- ggplot() + 
  stat_function(fun = mef_exp_y, col = "dodgerblue3", size = 0.7, xlim = c(0, 15)) +
  theme_bw() + scale_y_continuous(limit = c(9, 17), breaks = seq(9, 17, 2)) +
  labs(x = "time", y = "mean residual life") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15))

g_mef <- gridExtra::grid.arrange(g_mef_x, g_mef_y, ncol = 2)

#ggsave(plot = g_mef, "./Application/Fig3.pdf", height = 6, width = 15)
ggsave(plot = g_mef, "./Application/Fig3.eps", height = 6, width = 15)

mef = data.frame(time = c(1, 5,10,15,20,25,30),
                 mefw = mef_ray_x(c(1, 5,10,15,20,25,30)),
                 mefm = mef_exp_y(c(1, 5,10,15,20,25,30)))

xtable::xtable(t(mef))
