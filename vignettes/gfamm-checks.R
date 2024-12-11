

## ----load_libraries, message = FALSE-------------------------------------------------------------
library(GCFPCA)
library(tidyverse)
library(lme4)
library(refund)
library(mvtnorm)
library(splines)
library(tidyfun)
theme_set(theme_minimal())

plot_pffr <- function(m) {
  layout(t(1:3))
  plot(m, select = 1, shift = pffr_count$coefficients[1])
  lines(seq(0, 1, length.out = K), f_0(seq(0, 1, length.out = K)), col = 2, lty = 2)
  plot(m, select = 2)
  lines(seq(0, 1, length.out = K), f_1(seq(0, 1, length.out = K)), col = 2, lty = 2)
  re <- matrix(predict(m, type = "terms")[[3]], ncol = K)
  tfd(re, arg = m$pffr$yind) |> plot(ylab = "REs")
}

plot_pffr_eta <- function(m, sim) {
  layout(t(1:2))
  data <- tibble(
    pred = predict(m, type = "link") |> tfd(arg = m$pffr$yind),
    true = sim$df_gcfpca |> select(id, index, eta) |> tfd()
  ) |> pivot_longer(c(pred, true))
  ggplot(data) + geom_spaghetti(aes(y=value, col = name)) +
    facet_wrap(~name)
}

## ----pois_mod, message=FALSE, warning=FALSE------------------------------------------------------
# simulate data
library(splines)

I = 50
K = 1000

# simulate data and save results
set.seed(1133)

#set true fixed effects
theta_0 = rnorm(10 + 4, sd=1)
theta_1 = rnorm(10 + 4, sd=1)

#set true fixed effects
f_0 = function(s) bs(s, knots=seq(0.1, 0.9, len=10), Boundary.knots=c(0,1), intercept=TRUE) %*% theta_0
f_1 = function(s) bs(s, knots=seq(0.1, 0.9, len=10), Boundary.knots=c(0,1), intercept=TRUE) %*% theta_1

count_sim <- gcfpca_simu(I = I, K = K, family = "poisson",
                         beta0_true = f_0, beta1_true = f_1,
                         fe_case = 2, re_case = 2)

# use binwidth for 30 intervals to (roughly) compare with pffr model of similar
# coef-dimension using 25-40 basis functions:
gc_fpca_count_time <- system.time(
  count_model <- gc_fpca(formula = Y ~ X + (1|id),
                                   data = count_sim$df_gcfpca,
                                   binwidth = K/30,
                                   family = "poisson",
                                   pve = 0.95, npc = 4, periodicity = FALSE)
)
gc_fpca_count_time


count_sim_pffr <- tf_nest(count_sim$df_gcfpca, Y, eta, .id = id, .arg = index)
count_sim_pffr$Y_mat <- as.matrix(count_sim_pffr$Y)

pffr_count_re_time <- system.time({
  pffr_count_re <- pffr(Y_mat ~ X +
                          s(id, bs = "re", bs.yindex = list(k = 10)), # GC-FPCA uses 4 FPCs for this
                        data = count_sim_pffr, family = "poisson",
                           bs.int = list(k = 40), bs.yindex = list(k = 20),
                           yind = tf_arg(count_sim_pffr$Y),
                           algorithm = "bam", method= "fREML")
})
pffr_count_re_time

plot_pffr(pffr_count_re)
plot_pffr_eta(pffr_count_re, count_sim)

pffr_count_pcre_time <- system.time({
  pffr_count_pilot <- pffr(Y_mat ~ X, data = count_sim_pffr, family = "poisson",
                           bs.int = list(k = 40), bs.yindex = list(k = 20),
                           yind = tf_arg(count_sim_pffr$Y),
                           algorithm = "bam", method= "fREML")
  count_fpc_e <- resid(pffr_count_pilot) |> fpca.face(pve = .95)
  count_efuns <- count_fpc_e$efunctions
  count_evalues <- count_fpc_e$evalues
  pffr_count_pcre <- pffr(Y_mat ~ X +
                       pcre(id = id, efunctions = count_efuns,
                            evalues = count_evalues, yind = tf_arg(count_sim_pffr$Y)),
                     data = count_sim_pffr, family = "poisson",
                     bs.int = list(k = 40), bs.yindex = list(k = 20),
                     yind = tf_arg(count_sim_pffr$Y),
                     algorithm = "bam", method= "fREML", discrete = TRUE)
})
pffr_count_pcre_time

plot_pffr(pffr_count_pcre)
plot_pffr_eta(pffr_count_pcre, count_sim)


  ## ----plot_pFPCA_fe, echo = TRUE, fig.show='hold'-------------------------------------------------
  plot_df <- cbind.data.frame(sind = rep(seq(0, 1, length.out = K), 2),
                              betahat = c(data.matrix(count_model$betaHat)),
                              betatrue = c(f_0(seq(0, 1, length.out = K)),
                                           f_1(seq(0, 1, length.out = K))),
                              X = c(rep("Intercept", K),
                                    rep("X", K)),
                              CI_L_pw = c(data.matrix(count_model$CI_L_pw)),
                              CI_U_pw = c(data.matrix(count_model$CI_U_pw)),
                              CI_L_joint = c(data.matrix(count_model$CI_L_joint)),
                              CI_U_joint = c(data.matrix(count_model$CI_U_joint))) %>%
    mutate(X = factor(X, levels = c("Intercept",
                                    "X")))

  plot_df %>%
    ggplot(aes(x = sind, y = betahat)) +
    geom_ribbon(aes(ymin = CI_L_joint, ymax = CI_U_joint, fill = "CI Joint"), alpha = 0.5) +
    geom_ribbon(aes(ymin = CI_L_pw, ymax = CI_U_pw, fill = "CI"), alpha = 0.5) +
    geom_line(aes(color = "GCFPCA")) +
    geom_line(aes(x = sind, y = betatrue, color = "truth")) +
    scale_fill_manual(values = c("CI" = "black", "CI Joint" = "lightgray"), name = "Confidence Interval") +
    scale_color_manual(values = c("GCFPCA" = "darkblue", "truth" = "red"), name = "Confidence Interval") +
    # Adding a horizontal dotted line at y = 0
    geom_hline(yintercept = 0, linetype = "dotted") +
    # Setting x-axis labels to show time
    scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
    # Facet the plot by variable X, with 3 columns
    facet_wrap(~X, ncol = 3, scales = "free_y") +
    # Adding axis labels and title
    labs(x = "Functional Domain", y = "", fill = "Confidence Interval")


  ## ----plot_pFPCA_re, echo = TRUE, fig.show='hold'-------------------------------------------------
  flip_efunctions = function(data){
    if(data[["logit_mod"]] %*% data[["truth"]] < 0){
      data[["logit_mod"]] = -1 * data[["logit_mod"]]
    }
    data
  }

  phi_df = tibble(s = rep(seq(0, 1, length.out  = K), 4),
                  l = rep(paste0("eigenfunction ",1:4), each = K),
                  truth = c(count_sim$phi),
                  logit_mod = c(data.matrix(count_model$efunctions[, -1]))) %>%
    nest_by(l)

  new_phi = map(phi_df$data, flip_efunctions)
  phi_df$new_phi <- new_phi

  phi_df %>%
    select(-data) %>%
    unnest(new_phi) %>%
    ungroup() %>%
    mutate(logit_mod = logit_mod * sqrt(K)) %>%
    pivot_longer(truth:logit_mod, names_to = "model", values_to = "value") %>%
    ggplot(aes(s, value, group = model, color = model, linetype = model)) +
    geom_line() +
    facet_wrap(~l)


  ## ----plot_pFPCA_eta------------------------------------------------------------------------------
  count_sim$df_gcfpca %>%
    mutate(eta_hat = as.vector(count_model$etas)) %>%
    filter(id %in% c(1, 2, 8)) %>%
    ggplot(aes(index, eta)) +
    geom_line() +
    geom_line(aes(y = eta_hat), linetype = 2, color = "salmon") +
    facet_wrap(~id)



## ----df_long-------------------------------------------------------------------------------------
df_long <- readRDS(here::here("data/df_long.rds"))


## ----nhanes_gcfpca, warning=FALSE, message=FALSE-------------------------------------------------
gcfpca_start_t = Sys.time()
nhanes_gcfpca <- gc_fpca(formula = mims ~ age + gender +(1|id),
                         data = df_long,
                         binwidth = 100,
                         family = "binomial",
                         pve = 0.95, npc = 4, periodicity = TRUE)
gcfpca_end_t = Sys.time()
gcfpca_time_diff = as.double(difftime(gcfpca_end_t, gcfpca_start_t, units="mins"))
gcfpca_time_diff


## ----nhanes_pffr, warning=FALSE, message=FALSE---------------------------------------------------
I = length(unique(df_long$id))
K = length(unique(df_long$index))

df_pffr <- data.frame(Y = I(matrix(df_long$mims, I, K, byrow=TRUE)),
                      id = factor(unique(df_long$id)),
                      age = df_long$age[!duplicated(df_long$id)],
                      gender = df_long$gender[!duplicated(df_long$id)])

pffr_start_t = Sys.time()
pffr_model_re <- pffr(Y ~ age + gender + s(id, bs="re", bs.yindex = list(bs = "cc", k = 10)),
                   algorithm="bam", method="fREML", discrete=TRUE,
                   bs.yindex=list(bs = "cc", k=20),
                   bs.int=list(bs="cc", k=40),
                   data=df_pffr,
                   family="binomial", yind=1:1440)
pffr_end_t = Sys.time()
pffr_re_time_diff = as.double(difftime(pffr_end_t, pffr_start_t, units="mins"))
pffr_re_time_diff

pffr_start_t = Sys.time()
pffr_model_pilot <- pffr(Y ~ age + gender,
                      algorithm="bam", method="fREML", discrete=TRUE,
                      bs.yindex=list(bs = "cc", k = 20),
                      bs.int=list(bs="cc", k = 40),
                      data=df_pffr,
                      family="binomial", yind=1:1440)
nhanes_fpc_e <- resid(pffr_model_pilot) |> fpca.face(pve = .9) #need 11 FPCs for 90% "PVE"
nhanes_efuns <- nhanes_fpc_e$efunctions
nhanes_evalues <- nhanes_fpc_e$evalues

system.time(pffr_model <- pffr(Y ~ age + gender +
                                 pcre(id = id, efunctions = nhanes_efuns,
                                      evalues = nhanes_evalues, yind = 1:1440),
                               algorithm="bam", method="fREML", discrete = TRUE,
                               bs.yindex=list(bs = "cc", k=30),
                               bs.int=list(bs="cc", k = 30),
                               data=df_pffr,
                               family="binomial", yind=1:1440))
pffr_end_t = Sys.time()
(pffr_pc_time_diff = as.double(difftime(pffr_end_t, pffr_start_t, units="mins")))

sind = 1:1440
df_pred <- data.frame("yindex.vec" = sind,
                      id = as.numeric(levels(df_long$id))[1],
                      id.vec = as.numeric(levels(df_long$id))[1],
                      age = 1, gender = 1)
betahat_pffr <- mgcv::predict.gam(pffr_model, newdata=df_pred, type='iterms', se.fit=TRUE)
CI_L_pffr <- betahat_pffr[[1]] - 1.96 * betahat_pffr[[2]]
CI_U_pffr <- betahat_pffr[[1]] + 1.96 * betahat_pffr[[2]]

eta_pffr <- matrix(pffr_model$linear.predictors, nrow = 50, ncol = 1440)


## ----plot_nhanes_fe, echo = TRUE, fig.show='hold'------------------------------------------------
plot_df <- cbind.data.frame(sind = rep(seq(0, 1, length.out = K), 3),
                            betahat = c(data.matrix(nhanes_gcfpca$betaHat)),
                            betahat_pffr = c(data.matrix(betahat_pffr[[1]][, -4])),
                            X = c(rep("Intercept", K),
                                  rep("Age", K),
                                  rep("Gender", K)),
                            CI_L_pw = c(data.matrix(nhanes_gcfpca$CI_L_pw)),
                            CI_U_pw = c(data.matrix(nhanes_gcfpca$CI_U_pw)),
                            CI_L_joint = c(data.matrix(nhanes_gcfpca$CI_L_joint)),
                            CI_U_joint = c(data.matrix(nhanes_gcfpca$CI_U_joint)),
                            CI_L_pffr = c(CI_L_pffr[, -4]),
                            CI_U_pffr = c(CI_U_pffr[, -4])) %>%
  mutate(X = factor(X, levels = c("Intercept",
                                  "Age",
                                  "Gender")))

plot_df %>%
  ggplot(aes(x = sind, y = betahat)) +
  #geom_ribbon(aes(ymin = CI_L_joint, ymax = CI_U_joint, fill = "CI Joint"), alpha = 0.5) +
  geom_ribbon(aes(ymin = CI_L_pw, ymax = CI_U_pw, fill = "CI"), alpha = 0.3) +
  geom_ribbon(aes(ymin = CI_L_pffr, ymax = CI_U_pffr, fill = "CI pffr"), alpha = 0.5) +
  geom_line(aes(color = "GCFPCA")) +
  geom_line(aes(x = sind, y = betahat_pffr, color = "GFAMM")) +
  scale_fill_manual(values = c("CI" = "blue", "CI Joint" = "lightblue", "CI pffr" = "pink"), name = "Confidence Interval") +
  scale_color_manual(values = c("GCFPCA" = "blue", "GFAMM" = "red"), name = "Confidence Interval") +
  # Adding a horizontal dotted line at y = 0
  geom_hline(yintercept = 0, linetype = "dotted") +
  # Setting x-axis labels to show time
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  # Facet the plot by variable X, with 3 columns
  facet_wrap(~X, ncol = 3, scales = "free_y") +
  # Adding axis labels and title
  labs(x = "Functional Domain", y = "", fill = "Confidence Interval")


## ----plot_nhanes_eta-----------------------------------------------------------------------------
df_long %>%
  mutate(eta_hat = as.vector(nhanes_gcfpca$etas)) %>%
  mutate(eta_pffr = as.vector(eta_pffr)) %>%
  filter(id %in% c(63529, 67427, 82257)) %>%
  ggplot(aes(index, eta_pffr)) +
  geom_line() +
  geom_line(aes(y = eta_hat), linetype = 2, color = "salmon") +
  facet_wrap(~id)


## ----plot_nhanes_re, echo = TRUE, fig.show='hold'------------------------------------------------
phi_df = tibble(s = rep(seq(0, 1, length.out  = 1440), 4),
                l = rep(paste0("eigenfunction ", 1:4), each = 1440),
                logit_mod = c(data.matrix(nhanes_gcfpca$efunctions[, -1]))) %>%
  nest_by(l)

phi_df %>%
  unnest(data) %>%
  ungroup() %>%
  mutate(logit_mod = logit_mod * sqrt(1440)) %>%
  pivot_longer(logit_mod, names_to = "model", values_to = "value") %>%
  ggplot(aes(s, value, group = model, color = model, linetype = model)) +
  geom_line(color = "black") +
  facet_wrap(~l)

