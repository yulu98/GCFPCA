#' Simulate longitudinal functional response Y (matrix) and a scalar predictor X (vector)
#'
#' @import mvtnorm
#' @param family distribution of longitudinal functional data, including "gaussian", "binomial", "poisson".
#' @param I number of subjects.
#' @param K number of grid points on the functional domain.
#' @param beta0_true true intercept function.
#' @param beta1_true true fixed effect function.
#' @param fe_case Takes on values of 1 or 2. If case = 1, true fixed effect
#' coefficient is generated from a Bernoulli(0.5) distribution .
#' If case = 2 true fixed effect, coefficient is from a N(0,1) distribution.
#' @param re_case Takes on values of 1 or 2. If case = 1 true eigenfunctions
#' are based on alternating sine and cosines. If case = 2 true eigenfunctions are
#' based on sqrt functions.
#' @export
#' @return a data frame containing generated predictors and longitudinal functional outcomes.

gcfpca_simu <- function(family = "gaussian", I = 100, K = 100,
                        beta0_true, beta1_true,
                        fe_case = 1, re_case = 1){

  if(!(family %in% c("binomial", "poisson", "gaussian"))){
    stop('family must be "binomial", "poisson", or "gaussian" for simulated data.')
  }

  # functional domain, equally spaced grid on 0,1
  index <- seq(0, 1, len=K)

  ## fixed effect
  # generate fixed effect covariates
  if(fe_case == 1){
    X_i <- rbinom(I, size=1, prob=0.5)
  }else{
    X_i <- rnorm(I)
  }


  ## random effects
  # 4 FPCs used in simulating data
  L = 4
  # eigenfunctions \phi_k(s) evaluated on the observed grid
  if(re_case == 1){
    phi <- sqrt(2) * cbind(sin(2*pi*index), cos(2*pi*index),
                           sin(4*pi*index), cos(4*pi*index))
  }else if(re_case == 2){
    phi <- cbind(rep(1, K), sqrt(3)*(2*index-1),
                 sqrt(5)*(6*index^2-6*index+1),
                 sqrt(7)*(20*index^3-30*index^2+12*index-1))
  }
  # eigenvalues \lambda
  lambda <- 0.5^(0:(L-1))
  # subject-specific weights/coefficients \xi_ik
  xi <- matrix(rnorm(I*4), I, L)
  xi <- xi %*% diag(sqrt(lambda))
  b_i <- xi %*% t(phi); # subject-specific random effect of size I by K

  ## linear predictor \eta_i(s)
  eta_i <- t(vapply(1:I, function(x){
    beta0_true(index) + beta1_true(index)*X_i[x] + b_i[x,]
  }, numeric(K)))

  ## response Y_i(s)
  if(family == "binomial"){
    Y_i <- matrix(rbinom(I*K, size=1, prob=plogis(eta_i)),
                        I, K, byrow=FALSE)
  }else if(family == "poisson"){
    Y_i <- matrix(rpois(I*K, lambda=exp(eta_i)),
                  I, K, byrow=FALSE)
  }else if(family == "gaussian"){
    Y_i <-
    Y <- X + sigma*matrix(rnorm(N*J),N,J)
  }

  ## create data matrix in long format for estimation
  df <- data.frame(id = as.factor(rep(1:I, each = K)),
                   index = rep(index, I),
                   Y = as.vector(t(Y_i)),
                   eta = as.vector(t(eta_i)),
                   X = rep(X_i, each = K))

  list(
    df_gcfpca = df,
    phi = phi,
    lambda = lambda,
    scores = xi
  )
}

