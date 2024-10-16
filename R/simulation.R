library(splines)

I = 100
K = 100

# simulate data and save results
set.seed(1134)

#set true fixed effects
theta_0 = rnorm(10 + 4, sd=1)
theta_1 = rnorm(10 + 4, sd=1)

#set true fixed effects
f_0 = function(s) bs(s, knots=seq(0.1, 0.9, len=10), Boundary.knots=c(0,1), intercept=TRUE) %*% theta_0
f_1 = function(s) bs(s, knots=seq(0.1, 0.9, len=10), Boundary.knots=c(0,1), intercept=TRUE) %*% theta_1

simu_data <- gcfpca_simu(I = 100, K = 100, family = "poisson",
                         beta0_true = f_0, beta1_true = f_1)

gcfpca_fit <- gc_fpca(formula = Y ~ X + (1|id),
                      data = simu_data$df_gcfpca,
                      binwidth = 10,
                      family = "poisson",
                      pve = 0.95, npc = 4)

plot_df <- cbind.data.frame(sind = rep(seq(0, 1, length.out = K), 2),
                            betahat = c(data.matrix(gcfpca_fit$betaHat)),
                            X = c(rep("Intercept", K),
                                  rep("X", K)),
                            CI_L_pw = c(data.matrix(gcfpca_fit$CI_L_pw)),
                            CI_U_pw = c(data.matrix(gcfpca_fit$CI_U_pw)),
                            CI_L_joint = c(data.matrix(gcfpca_fit$CI_L_joint)),
                            CI_U_joint = c(data.matrix(gcfpca_fit$CI_U_joint))) %>%
  mutate(X = factor(X, levels = c("Intercept",
                                  "X")))

plot_df %>%
  ggplot(aes(x = sind, y = betahat)) +
  geom_ribbon(aes(ymin = CI_L_joint, ymax = CI_U_joint, fill = "CI Joint"), alpha = 0.5) +
  geom_ribbon(aes(ymin = CI_L_pw, ymax = CI_U_pw, fill = "CI"), alpha = 0.5) +
  geom_line(color = "darkblue") +
  scale_fill_manual(values = c("CI" = "black", "CI Joint" = "lightgray"), name = "Confidence Interval") +
  # Adjusting the theme for a clean plot appearance
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(),

    # Customize axis text and title
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),

    # Customize strip text for facets
    strip.text = element_text(size = 20),
    strip.background = element_rect(colour = "white", fill = "white"),
    # Legend customization
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15),
    legend.key = element_rect(fill = "white"),
    # Other elements
    legend.title = element_text(size = 15)
  ) +
  # Adding a horizontal dotted line at y = 0
  geom_hline(yintercept = 0, linetype = "dotted") +
  # Setting x-axis labels to show time
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  # Facet the plot by variable X, with 3 columns
  facet_wrap(~X, ncol = 3, scales = "free_y") +
  # Adding axis labels and title
  labs(x = "Functional Domain", y = "", fill = "Confidence Interval")


ef_plot <- gcfpca_fit$efunctions
ef_plot %>%
  pivot_longer(phi1:phi4, names_to = "model", values_to = "value") %>%
  mutate(model = factor(model, levels = paste0("phi", 1:4))) %>%
  ggplot(aes(x = index, y = value)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_wrap(~ model, ncol = 4, scales = "free_y") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(legend.title = element_blank(),
        legend.position = c(0.9, -0.17), # Adjust position to top-right of the left plot
        legend.direction = "horizontal",
        legend.text = element_text(size = 15),
        legend.key = element_rect(fill = "white")) +
  theme(strip.text = element_text(size = 20),
        strip.background = element_rect(fill=NA),
        panel.border = element_rect(fill = NA, colour = "black")) +
  theme(panel.spacing = unit(0.5, "lines"), # Adjust space between panels if needed
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.line.y.left = element_line(colour = "black")) + # Add border to each panel
  labs(x = "Time of Day",
       y = "")
