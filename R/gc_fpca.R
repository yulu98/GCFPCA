#' Generalized Conditional Functional Principal Components Analysis (GC-FPCA)
#'
#' This is the main function for the GCFPCA package.
#'
#' @import dplyr
#' @importFrom stats coef predict binomial lm median as.formula
#' @importFrom refund fpca.face
#' @importFrom lme4 glmer
#' @importFrom utils txtProgressBar setTxtProgressBar data
#' @importFrom mgcv bam predict.bam
#' @importFrom mvtnorm rmvnorm
#' @importFrom magrittr %>%
#'
#' @param formula A formula specifying the GLMM model.
#' @param data A data frame with columns: `id` (subject), `index` (time or index), and `Y` (outcome variable).
#' @param binwidth Numeric; controls the width of the bins for local fitting. Defaults to 10.
#' @param pve Numeric; proportion of variance explained, used to choose the number of principal components unless `npc` is specified.
#' @param npc Numeric; number of smooth PCs to extract. If `NULL`, `npc` is chosen based on `pve`.
#' @param family A family object for specifying the distribution and link function for the GLMM. Defaults to `binomial`.
#' @param ... Additional arguments passed to the `mgcv::bam` function.
#'
#'
#' @export
#' @return An object of class `gc_fpca` containing:
#' \item{betaHat}{Estimated fixed effects coefficient functions}
#' \item{betaHat.sd}{Estimated fixed effects coefficient standard deviation along the functional domain}
#' \item{efunctions}{\eqn{K \times npc} matrix of estimated FPC basis functions.}
#' \item{scores}{\eqn{I \times npc} matrix of estimated FPC scores.}
#' \item{mu}{Estimated population-level mean function.}
#' \item{Yhat}{FPC approximation of subject-specific means.}
#' \item{family}{The family used for GLMM.}
#'
#' @examples
#' # Simulate data and fit GC-FPCA
#' df_list <- gen_data(N = 100, J = 50, run_num = 1)
#' gcfpca_mod <- gc_fpca(df_list, bin_width = 10, family = binomial())
gc_fpca <- function(formula, data, binwidth = 10, family = "gaussian",
                    pve = NULL, npc = NULL, periodicity = FALSE,
                    ...){
  ## Organize the input
  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)
  Y_name <- all.vars(formula)[1]  # Extracts 'Y'
  X_name <- all.vars(formula)[-c(1, length(all.vars(formula)))] # Extracts 'X'
  p = length(X_name)
  I = length(unique(data$id))
  K = length(unique(data$index))
  index = sort(unique(data$index))
  if(K/binwidth < 10){
    binwidth = K/10
    message(paste0("binwidth should be no more than K/10. Converting to a new binwidth of ", binwidth, "."))
  }

  ##########################################################################################
  ## Step 1, 2
  ##########################################################################################

  fit_ls  <- vector(mode = "list", length = K)
  for(k in 1:K){
    # Step 1: Create data bins.
    if(periodicity == TRUE){
      ind_k <- (k - floor(binwidth/2)):(k + floor(binwidth/2)) %% K
      ind_k[ind_k == 0] <- K
      data_k <- data %>% filter(index %in% index[ind_k])
    }else{
      ind_k <- (k - floor(binwidth/2)):(k + floor(binwidth/2))
      ind_k <- ind_k[ind_k > 0]
      data_k <- data %>% filter(index %in% index[ind_k])
    }


    # Step 2: Fit local GLMMs; extract linear predictor estimates.
    fit_k <- glmer(
      formula,
      data = data_k, family = family,
      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
    )

    fit_ls[[k]] <- data.frame("id" = 1:I,
                              "eta_i" = coef(fit_k)$id[[1]]-fixef(fit_k)[1],
                              "index_bin" = k)
  }

  ##########################################################################################
  ## Step 3
  ##########################################################################################
  # Combine list of fitted data into a single dataframe and arrange by id and time
  fit_df <- fit_ls %>%
    bind_rows() %>%
    arrange(id, index_bin)

  # Pivot the dataframe to a wide format for FPCA input, with `eta_i` values spread across time points
  fit_df_wide <- fit_df %>%
    select(id, index_bin, eta_i) %>%
    pivot_wider(id_cols = id, names_from = index_bin, values_from = eta_i, names_prefix = "S")

  # Apply FPCA on the wide matrix of eta_i values
  if(K/binwidth < 20){
    knots <- ceiling(K/binwidth)
  }else{
    knots <- 20
  }

  fpca_latent <- fpca.face(as.matrix(select(fit_df_wide, -id)),
                           pve = pve,
                           npc = npc,
                           argvals =1:K,
                           knots = knots,
                           lower=0)

  # Extract eigenfunctions and associate them with corresponding index (time point index)
  efun <- data.frame(fpca_latent$efunctions) %>%
    mutate(index = index)  # Ensure `index` aligns with the unique time indices

  # Merge the eigenfunctions back into the original dataframe by `index`
  df_w_ef <- data %>%
    left_join(efun, by = "index")

  ##########################################################################################
  ## Step 4
  ##########################################################################################
  # Dynamically change column names based on number of eigenfunctions
  df_w_ef <- df_w_ef %>%
    rename_with(~ gsub("\\.", "phi1", .)) %>%
    rename_with(~ gsub("X(?=[0-9])", "phi", ., perl = TRUE))

  # Count how many columns start with 'phi', representing the number of eigenfunctions
  phi_num <- sum(startsWith(colnames(df_w_ef), "phi"))

  # Dynamically construct GAM formula based on number of eigenfunctions
  # If there are no eigenfunctions (phi_num = 0), exclude them from the formula
  add_ef <- case_when(
    phi_num == 0 ~ "",
    phi_num != 0 ~ paste0(
      # Create a list of GAM terms for each eigenfunction and collapse into a string
      lapply(1:phi_num, function(x) paste0("s(id, by = phi", x, ", bs = 're')")),
      collapse = " + "
    )
  )

  # Build the final dynamic formula by incorporating the eigenfunction terms
  dynamic_form <- as.formula(paste0(
    Y_name, " ~ s(index, k = ", 30, ", bs = 'cr') + ",
    paste0(paste0("s(index, by = ", X_name, ", k = ", 30, ", bs = 'cr') + "), collapse = ""),
    add_ef
  ))

  # Fit the GAM model using the dynamic formula
  fit_gcfpca <- mgcv::bam(
    formula = dynamic_form,
    method = "fREML",
    data = df_w_ef,
    discrete = TRUE,
    family = family,
    gc.level = 1
  )

  ###################################################################
  # Organize the results
  ###################################################################
  ############ Fixed effects with pointwise and joint CI ############
  sample_id = as.numeric(levels(data$id))[1]
  df_pred = paste0("data.frame(index = index, ",
                   paste0(paste0("phi", 1:phi_num), collapse = " = 0, "),
                   " = 0, id = ", sample_id, ", ",
                   paste0(X_name, collapse = " = 1, "), "= 1)")
  df_pred = eval(parse(text = df_pred))

  pred_fit = predict(fit_gcfpca, newdata=df_pred, se.fit=TRUE, type = "terms")

  Vb <- vcov(fit_gcfpca)
  X_pred <- predict(fit_gcfpca, newdata=df_pred, type="lpmatrix")

  s_MIN_ind <- grep("s\\(index\\)", colnames(X_pred))
  intercept_ind <- which(colnames(X_pred) == "(Intercept)")
  s_MIN_ind <- c(intercept_ind, s_MIN_ind)
  X_pred_s_MIN <- X_pred[, s_MIN_ind]
  cov_s_MIN <- X_pred_s_MIN %*% Vb[s_MIN_ind, s_MIN_ind] %*% t(X_pred_s_MIN)

  cov_s_X_MIN <- list()
  for(x in X_name){
    s_min_X_ind <- grep(paste0("s\\(index\\):", x), colnames(X_pred))
    X_pred_s_X_MIN <- X_pred[, s_MIN_ind]
    cov_s_X_MIN[[paste0(x)]] <- X_pred_s_MIN %*% Vb[s_MIN_ind, s_MIN_ind] %*% t(X_pred_s_MIN)
  }

  pred_val = pred_fit$fit[, 1:(p+1)]
  pred_val[, 1] = pred_val[, 1] + fit_gcfpca$coefficients[1]
  pred_sd = pred_fit$se.fit[, 1:(p+1)]
  pred_sd[, 1] = sqrt(diag(cov_s_MIN))

  CI_L_pw = pred_val - 1.96 * pred_sd
  CI_U_pw = pred_val + 1.96 * pred_sd

  ## obtain qn to construct joint CI
  qn <- rep(0, length = p+1)
  sigma <- c(list(intercept = cov_s_MIN), cov_s_X_MIN)
  N <- 10000 ## sample size in simulation-based approach
  for(i in 1:(p+1)){
    x_sample <- rmvnorm(N, mean = pred_val[, i], sigma = sigma[[i]])
    un <- rep(NA, N)
    for(j in 1:N){
      un[j] <- max(abs((x_sample[j,] - pred_val[, i])/pred_sd[, i]))
    }
    qn[i] <- quantile(un, 0.95)
  }

  CI_L_joint = pred_val -  sweep(pred_sd, 2, qn, `*`)
  CI_U_joint = pred_val + sweep(pred_sd, 2, qn, `*`)

  fixed_effect <- list(betaHat = pred_val,
                       betaHat.sd = pred_sd,
                       CI_L_pw = CI_L_pw,
                       CI_U_pw = CI_U_pw,
                       CI_L_joint = CI_L_joint,
                       CI_U_joint = CI_U_joint)
  ############# random effects ############
  # report the random effect in the format of FACE results
  efunctions <- df_w_ef %>% filter(id == sample_id) %>% select(index, starts_with("phi"))
  score_hat <- coef(fit_gcfpca)
  scores <- matrix(score_hat[grep("s\\(id\\):phi",names(score_hat))], I, npc)

  ################# eta ###################
  etas <- matrix(predict(fit_gcfpca, type = "link"),
                 nrow = I, ncol = K)

  return(c(fixed_effect,
           list(efunctions = efunctions,
                scores = scores,
                etas = etas),
           model = list(fit_gcfpca))) # !
}
