library(doParallel)
library(foreach)
library(here)
library(splines)
library(GCFPCA)
library(refund)
library(tidyfun)
library(tidyverse)
library(patchwork)
library(lme4)


plot_pffr <- function(m) {
  layout(t(1:3))
  plot(m, select = 1,  scheme = 1)
  lines(seq(0, 1, length.out = K), f_0(seq(0, 1, length.out = K)), col = 2)
  plot(m, select = 2, scheme = 1)
  lines(seq(0, 1, length.out = K), f_1(seq(0, 1, length.out = K)), col = 2)
  re <- matrix(predict(m, type = "terms")[[3]], ncol = K)
  tfd(re) |> plot(ylab = "REs")
}

plot_pffr_eta <- function(m, sim) {
  data <- tibble(
    true = sim$df_gcfpca |> select(id, index, eta) |> tfd(),
    pred = predict(m, type = "link") |> tfd(arg = tf_arg(true))
  ) |>
    pivot_longer(c(pred, true))
  ggplot(data) +
    geom_spaghetti(aes(y = value, col = name)) +
    facet_wrap(~name)
}

plot_gcfpca <- function(count_model) {
  plot_df <- cbind.data.frame(
    sind = rep(seq(0, 1, length.out = K), 2),
    betahat = c(data.matrix(count_model$betaHat)),
    betatrue = c(
      f_0(seq(0, 1, length.out = K)),
      f_1(seq(0, 1, length.out = K))
    ),
    X = c(
      rep("Intercept", K),
      rep("X", K)
    ),
    CI_L_pw = c(data.matrix(count_model$CI_L_pw)),
    CI_U_pw = c(data.matrix(count_model$CI_U_pw)),
    CI_L_joint = c(data.matrix(count_model$CI_L_joint)),
    CI_U_joint = c(data.matrix(count_model$CI_U_joint))
  ) %>%
    mutate(X = factor(X, levels = c(
      "Intercept",
      "X"
    )))

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
}

plot_gcfpca_eta <- function(m, sim) {
  sim$df_gcfpca %>%
    mutate(eta_hat = as.vector(m$etas)) %>%
    pivot_longer(cols = c(eta, eta_hat)) |>
    # filter(id %in% c(1, 2, 8)) %>%
    ggplot(aes(x = index, y = value, group = id, col = name)) +
    geom_line() +
    facet_wrap(~ rev(name))
}


# register parallel backend
num_cores <- 5
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# define parameter grids
I_vals <- c(10, 20, 40, 80, 160) # example values
K_vals <- c(100, 200, 400, 1000) # example values

foreach(I = I_vals, .export = ls(),
        .packages = c("splines","GCFPCA","refund", "tidyverse", "tidyfun", "patchwork","here", "lme4", "mvtnorm")) %:%
  foreach(K = K_vals, .export = ls(),
          .packages = c("splines","GCFPCA","refund", "tidyverse", "tidyfun", "patchwork","here", "lme4", "mvtnorm")) %dopar% {

    cat("starting I:", I, "K:", K, "\n")
    theme_set(theme_minimal())

    theta_0 <- rnorm(10 + 4, sd = 1)
    theta_1 <- rnorm(10 + 4, sd = 1)
    theta_0 <- scale(theta_0, scale = FALSE)
    theta_1 <- scale(theta_1, scale = FALSE)

    f_0 <- function(s) bs(s, knots = seq(0.1, 0.9, len = 10),
                          Boundary.knots = c(0, 1), intercept = TRUE) %*% theta_0
    f_1 <- function(s) bs(s, knots = seq(0.1, 0.9, len = 10),
                          Boundary.knots = c(0, 1), intercept = TRUE) %*% theta_1

    # simulate data
    count_sim <- gcfpca_simu(
      I = I, K = K, family = "poisson",
      beta0_true = f_0, beta1_true = f_1,
      fe_case = 2, re_case = 2
    )

    # fit gc_fpca model
    gc_fpca_count_time <- system.time(
      count_model <- gc_fpca(
        formula = Y ~ X + (1 | id),
        data = count_sim$df_gcfpca,
        binwidth = K / 30,
        family = "poisson",
        pve = 0.95, npc = 4, periodicity = FALSE
      )
    )

    # fit pffr with RE
    count_sim_pffr <- tf_nest(count_sim$df_gcfpca, Y, eta, .id = id, .arg = index)
    count_sim_pffr$Y_mat <- as.matrix(count_sim_pffr$Y)
    pffr_count_re_time <- system.time({
      pffr_count_re <- pffr(
        Y_mat ~ X +
          s(id, bs = "re", bs.yindex = list(k = 10, bs = "cr")),
        data = count_sim_pffr, family = "poisson",
        bs.int = list(k = 30, bs = "cr"), bs.yindex = list(k = 30, bs = "cr"),
        yind = tf_arg(count_sim_pffr$Y),
        algorithm = "bam", method = "fREML", discrete = TRUE,
      )
    })

    # fit pffr with pcre
    pffr_count_pcre_time <- system.time({
      pffr_count_pilot <- pffr(
        Y_mat ~ X,
        data = count_sim_pffr, family = "poisson",
        bs.int = list(k = 30, bs = "cr"), bs.yindex = list(k = 30, bs = "cr"),
        yind = tf_arg(count_sim_pffr$Y),
        algorithm = "bam", method = "fREML", discrete = TRUE,
      )
      count_fpc_e <- fpca.face(resid(pffr_count_pilot), npc = 4, pve = .95, knots = 20, lower = 0)
      count_efuns <- count_fpc_e$efunctions
      count_evalues <- count_fpc_e$evalues
      pffr_count_pcre <- pffr(
        Y_mat ~ X +
          pcre(id = id, efunctions =  count_efuns,
               evalues = count_evalues, yind = tf_arg(count_sim_pffr$Y)),
        data = count_sim_pffr, family = "poisson",
        bs.int = list(k = 30, bs = "cr"), bs.yindex = list(k = 30, bs = "cr"),
        yind = tf_arg(count_sim_pffr$Y),
        algorithm = "bam", method = "fREML", discrete = TRUE
      )
    })

    # save timing info
    results <- data_frame(
      method = c("gcfpca", "re", "pcre"),
      I = I,
      K = K,
      times = c(gc_fpca_count_time[3], pffr_count_re_time[3], pffr_count_pcre_time[3])
    ) |> bind_cols(AIC(count_model$model, pffr_count_re, pffr_count_pcre))

    saveRDS(results, file = here("attic", paste0("results_I", I, "_K", K, ".rds")))

    # produce and save plots
    pdf(file = here("attic", paste0("plots_I", I, "_K", K, ".pdf")))
    print(plot_gcfpca(count_model) / plot_gcfpca_eta(count_model, count_sim))
    print(plot_pffr(pffr_count_re))
    print(plot_pffr_eta(pffr_count_re, count_sim))
    print(plot_pffr(pffr_count_pcre))
    print(plot_pffr_eta(pffr_count_pcre, count_sim))
    dev.off()
    cat("DONE with I:", I, "K:", K, "\n")
    NULL
  }

stopCluster(cl)

#------------------------------

library(here)
library(stringr)
library(dplyr)

# list all fit_times files
results_files <- list.files(path = here("attic"), pattern = "^results_I\\d+_K\\d+\\.rds$", full.names = TRUE)

# load and combine into dataframe
results_df <- do.call(rbind, lapply(results_files, readRDS))

ggplot(results_df) +
  geom_line(aes(x = I, y = times, col = method)) +
  facet_grid(~ K, labeller = label_both)
ggsave(here("attic/sim_gfamm_check-times-K.pdf"))

ggplot(results_df) +
  geom_line(aes(x = K, y = times, col = method)) +
  facet_grid(~ I, labeller = label_both)
ggsave("attic/sim_gfamm_check-times-I.pdf")

ggplot(results_df) +
  geom_point(aes(y = factor(I):factor(K), x = AIC, col = method), position = position_dodge(width = .5)) +
  scale_x_log10()
ggsave("attic/sim_gfamm_check-AIC.pdf")
