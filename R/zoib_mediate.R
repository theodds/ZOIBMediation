#' Perform causal mediation analysis in the zero-one inflated beta regression
#' model
#'
#' @param zoib_fit Model fit with the bayes_zoib function
#' @param data Original data used to get zoib_fit
#' @param formula_m formula for the mediator used in constructing zoib_fit
#' @param formula_y formula for the outcome used in constructign zoib_fit
#' @param g_comp_thin Thinning interval used for the g-computation, defauls to 1
#' @param print_interval How many iterations before printing progress?
#' @param num_copy Number of pseudo-datasets to use when performing g-computation
#' @param med_name Name of the mediator in the original data
#' @param trt_name Name of the treatment in the original data
#' @param out_name Name of the outcome in the original data
#'
#' @return A data.frame containing samples from the posterior of the direct (zeta),
#' indirect (delta), and total effects (tau) for different values of the treatment,
#' as well as the expected value of the mediator for the different treatments.
#'
#' The output can be fed into zoib_sensitivity and plot_zoib_sensitivity, or
#' post-processed directly using the data.frame.
zoib_mediate <- function(zoib_fit,
                         data,
                         formula_m = ~ 1,
                         formula_y = ~ 1,
                         g_comp_thin = 1,
                         print_interval = 100,
                         num_copy = 2,
                         med_name = "M",
                         trt_name = "A",
                         out_name = "Y") {



  X_m <- model.matrix(formula_m, data)
  X_y <- model.matrix(formula_y, data)
  N   <- nrow(X_m)

  data_0 <- data; data_0[[trt_name]] <- 0
  data_1 <- data; data_1[[trt_name]] <- 1
  X_m_0 <- model.matrix(formula_m, data_0)
  X_m_1 <- model.matrix(formula_m, data_1)
  X_y_0 <- model.matrix(formula_y, data_0)
  X_y_1 <- model.matrix(formula_y, data_1)

  g_comp_iters <- seq(from = 1, to = length(as.matrix(zoib_fit, pars = "lp__")), by = g_comp_thin)

  # samples <- sampling(zoib_model, data = stan_data, pars = c("beta_mediator", "beta_outcome"))
  samples_mat <- as.matrix(zoib_fit, pars = c("beta_mediator","beta_outcome"))|>
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

  delta_0 <- numeric(length(g_comp_iters))
  delta_1 <- numeric(length(g_comp_iters))
  zeta_0  <- numeric(length(g_comp_iters))
  zeta_1  <- numeric(length(g_comp_iters))
  tau     <- numeric(length(g_comp_iters))
  for(r in 1:length(g_comp_iters)) {
    i         <- g_comp_iters[r]

    for(k in 1:num_copy) {
      ## Simulate the mediators

      U         <- runif(N)
      V         <- runif(N)
      fm        <- function(x) sim_zoib_copula_df(U, x, beta_m_alpha[i,],
                                                  beta_m_gamma[i,], beta_m_mu[i,],
                                                  beta_m_phi[i,])
      M_0       <- fm(X_m_0)
      M_1       <- fm(X_m_1)

      ## Compute the outcome means
      data_00 <- data_0; data_00[med_name] <- M_0
      data_01 <- data_0; data_01[med_name] <- M_1
      data_10 <- data_1; data_10[med_name] <- M_0
      data_11 <- data_1; data_11[med_name] <- M_1
      X_y_00 <- model.matrix(formula_y, data_00)
      X_y_01 <- model.matrix(formula_y, data_01)
      X_y_10 <- model.matrix(formula_y, data_10)
      X_y_11 <- model.matrix(formula_y, data_11)

      f      <- function(x) sim_zoib_copula_df(V, x, beta_y_alpha[i,],
                                               beta_y_gamma[i,],
                                               beta_y_mu[i,], beta_y_phi[i,])
      Y_00   <- f(X_y_00)
      Y_01   <- f(X_y_01)
      Y_10   <- f(X_y_10)
      Y_11   <- f(X_y_11)

      delta_0[r] <- sum(omega[i,] * (Y_01 - Y_00))
      delta_1[r] <- sum(omega[i,] * (Y_11 - Y_10))
      zeta_0[r]  <- sum(omega[i,] * (Y_10 - Y_00))
      zeta_1[r]  <- sum(omega[i,] * (Y_11 - Y_01))
      tau[r]     <- sum(omega[i,] * (Y_11 - Y_00))

      my_df <- tibble(
        Iteration = r,
        CopyID    = k,
        Param     = c(delta_0[r], delta_1[r], zeta_0[r], zeta_1[r], tau[r]),
        ParamName = c("delta_0", "delta_1", "zeta_0", "zeta_1", "tau")
      )

      out_dfs[[length(out_dfs) + 1]] <- my_df
    }

    if(r %% print_interval == 0) cat(paste("Finishing iteration", r, "\n"))

  }

  out <- do.call(rbind, out_dfs)
  rm(out_dfs)
  return(effect_df = out)

}



