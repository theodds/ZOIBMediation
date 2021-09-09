#' Sensitivity analysis for ZOIB mediation model
#'
#' Performs the sensitivity analysis outlined by Rene et al. (2021) Causal
#' Mediation and Sensitivity Analysis for Mixed-Scale Data - see this paper for
#' additional details. This sensitivity analysis shifts unobserved potential
#' outcomes Y(a, m) either on the linear or logit scale by a factor lambda *
#' (M(a) - m), where lambda is an unidentified sensitivity parameter.
#'
#' @param data Data used to fit the original model
#' @param samples Model fit obtained by bayes_zoib
#' @param formula_m Formula for the mediator when fitting the model
#' @param formula_y Formula for the outcome when fitting the model
#' @param g_comp_thin Do we want to sample on all iterations (thin = 1) or half (thin = 2), etc
#' @param print_interval Number of iterations before printing an update on progress
#' @param num_copy Number of pseudo-datasets when performing g-computation
#' @param med_name Name of the mediator
#' @param trt_name Name of the treatment
#' @param out_name Name of the outcome
#' @param rho Correlation in Gaussian copula (only relevant when doing sensitivity on the logit scale)
#' @param lambda_grid The values of lambda we want to look at
#' @param logit_scale Whether to operate on the logit scale or not.
#'
#' @return Samples of the mediation effects
#' (delta(0), delta(1), zeta(0), zeta(1), tau) for each value of lambda, stored
#' as a data.frame with the iteration, id of the simulated dataset, parameter,
#' and value of lambda given.
#'
#' This can be worked with directly for inference, but it is probably more efficient to
#' visualize the results with plot_zoib_sensitivity.
#' @export
zoib_sensitivity <- function(data, samples,
                             formula_m = ~ 1,
                             formula_y = ~ 1,
                             g_comp_thin = 1,
                             print_interval = 100,
                             num_copy = 2,
                             med_name,
                             trt_name,
                             out_name,
                             rho = 1,
                             lambda_grid,
                             logit_scale = TRUE
) {

  X_m <- model.matrix(formula_m, data)
  N   <- nrow(X_m)

  g_comp_iters <- floor(seq(from = 1, to = length(as.matrix(samples, pars = "lp__")),
                      by = g_comp_thin))

  data_0 <- data; data_0[[trt_name]] <- 0
  data_1 <- data; data_1[[trt_name]] <- 1
  X_m_0 <- model.matrix(formula_m, data_0)
  X_m_1 <- model.matrix(formula_m, data_1)
  lambda_mat <- matrix(rep(lambda_grid, N), nrow = N, byrow = TRUE)

  samples_mat <- as.matrix(samples, pars = c("beta_mediator","beta_outcome"))|>
    as.data.frame()

  beta_m_alpha <- samples_mat |> select(contains("mediator")) |>
    select(contains(",1]"))  |> as.matrix()
  beta_m_gamma <- samples_mat |> select(contains("mediator")) |>
    select(contains(",2]"))  |> as.matrix()
  beta_m_mu <- samples_mat |> select(contains("mediator")) |>
    select(contains(",3]"))  |> as.matrix()
  beta_m_phi <- samples_mat |> select(contains("mediator")) |>
    select(contains(",4]"))  |> as.matrix()

  beta_y_alpha <- samples_mat |> select(contains("outcome")) |>
    select(contains(",1]")) |> as.matrix()
  beta_y_gamma <- samples_mat |> select(contains("outcome")) |>
    select(contains(",2]"))  |> as.matrix()
  beta_y_mu <- samples_mat |> select(contains("outcome")) |>
    select(contains(",3]"))  |> as.matrix()
  beta_y_phi <- samples_mat |> select(contains("outcome")) |>
    select(contains(",4]"))  |> as.matrix()

  omega <- MCMCpack::rdirichlet(n = nrow(beta_m_alpha),
                                alpha = rep(1, nrow(X_m)))
  out_dfs <- list()

  Sigma <- rbind(c(1, rho), c(rho, 1))
  for(r in 1:length(g_comp_iters)) {
    i         <- g_comp_iters[r]

    for(k in 1:num_copy) {
      ## Simulate the mediators

      U <- pnorm(MASS::mvrnorm(n = N, mu = rep(0,2), Sigma = Sigma))

      fm        <- function(u, x) sim_zoib_copula_df(u, x, beta_m_alpha[i,],
                                                  beta_m_gamma[i,], beta_m_mu[i,],
                                                  beta_m_phi[i,])
      M_0       <- fm(U[,1], X_m_0)
      M_1       <- fm(U[,2], X_m_1)

      J         <- length(lambda_grid)
      M_0_mat   <- matrix(rep(M_0, J), nrow = N, ncol = J)
      M_1_mat   <- matrix(rep(M_1, J), nrow = N, ncol = J)

      ## Compute the outcome means
      data_00 <- data_0; data_00[med_name] <- M_0
      data_01 <- data_0; data_01[med_name] <- M_1
      data_10 <- data_1; data_10[med_name] <- M_0
      data_11 <- data_1; data_11[med_name] <- M_1
      X_y_00 <- model.matrix(formula_y, data_00)
      X_y_01 <- model.matrix(formula_y, data_01)
      X_y_10 <- model.matrix(formula_y, data_10)
      X_y_11 <- model.matrix(formula_y, data_11)

      f <- function(x) mean_zoib_df(X = x, beta_alpha = beta_y_alpha[i,],
                                    beta_gamma = beta_y_gamma[i,],
                                    beta_mu = beta_y_mu[i,],
                                    beta_phi = beta_y_phi[i,])

      Y_00   <- f(X_y_00) * matrix(1, nrow = N, ncol = length(lambda_grid))
      Y_01   <- f(X_y_01) * matrix(1, nrow = N, ncol = length(lambda_grid))
      Y_10   <- f(X_y_10) * matrix(1, nrow = N, ncol = length(lambda_grid))
      Y_11   <- f(X_y_11) * matrix(1, nrow = N, ncol = length(lambda_grid))

      if(logit_scale) {
        Y_10   <- expit(logit(Y_10) + lambda_mat * (M_1_mat - M_0_mat))
        Y_01   <- expit(logit(Y_01) + lambda_mat * (M_0_mat - M_1_mat))
      } else {
        Y_10 <- Y_10 + lambda_mat * (M_1_mat - M_0_mat)
        Y_01 <- Y_01 + lambda_mat * (M_0_mat - M_1_mat)
      }


      delta_0 <- colSums(omega[i,] * (Y_01 - Y_00))
      delta_1 <- colSums(omega[i,] * (Y_11 - Y_10))
      zeta_0  <- colSums(omega[i,] * (Y_10 - Y_00))
      zeta_1  <- colSums(omega[i,] * (Y_11 - Y_01))
      tau     <- colSums(omega[i,] * (Y_11 - Y_00))
      EM1     <- mean(M_1)
      EM0     <- mean(M_0)

      fdf        <- function(str, val)
        tibble(Iteration = r, CopyID = k, Param = val,
               ParamName = str, lambda = lambda_grid)
      delta_0_df <- fdf("delta_0", delta_0)
      delta_1_df <- fdf("delta_1", delta_1)
      zeta_0_df  <- fdf("zeta_0", zeta_0)
      zeta_1_df  <- fdf("zeta_1", zeta_1)
      tau_df     <- fdf("tau", tau)
      em1_df     <- fdf("EM1", EM1)
      em0_df     <- fdf("EM0", EM0)

      my_df <- rbind(delta_0_df, delta_1_df, zeta_0_df, zeta_1_df, tau_df, em1_df, em0_df)
      out_dfs[[length(out_dfs) + 1]] <- my_df
    }

    if(r %% print_interval == 0) cat(paste("Finishing iteration", r, "\n"))

  }

  out <- do.call(rbind, out_dfs)
  rm(out_dfs)
  return(out)
}
