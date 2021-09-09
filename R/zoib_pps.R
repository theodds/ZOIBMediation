#' Posterior predictions from ZOIB marginals
#'
#' Draws from the posterior predictive distribution of the mediator and outcome
#'
#' @param samples Model fit with the bayes_zoib function
#' @param thin Do we want to sample on all iterations (thin = 1) or half (thin = 2), etc.
#' @param data The data used to fit the model
#' @param formula_m Formula used to fit the mediator model
#' @param formula_y Formula used to fit the outcome model
#' @param med_name Name of the mediator
#' @param trt_name Name of the treatment
#' @param print_interval Number of iterations before printing an update on progress
#'
#' @return Samples from the posterior predictive distribution, stored as a data.frame
#' giving the subject id, iteration, outcome, and mediator.
zoib_ppc <- function(samples, thin, data, formula_m, formula_y,
                     med_name, trt_name, print_interval = 100) {

  iterations <- floor(seq(from = 1, to = length(as.matrix(samples, pars = "lp__")),
                    by = thin))

  X_m <- model.matrix(formula_m, data)
  N   <- nrow(X_m)

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

  out_samples <- list()
  for(r in 1:length(iterations)) {
    i  <- iterations[r]
    U  <- runif(N)
    V  <- runif(N)
    fm <- function(x) sim_zoib_copula_df(U, x, beta_m_alpha[i,],
                                         beta_m_gamma[i,], beta_m_mu[i,],
                                         beta_m_phi[i,])
    f  <- function(x) sim_zoib_copula_df(V, x, beta_y_alpha[i,],
                                         beta_y_gamma[i,],
                                         beta_y_mu[i,], beta_y_phi[i,])

    M  <- fm(X_m)
    data_new <- data; data_new[med_name] <- M
    X_y <- model.matrix(formula_y, data_new)
    Y <- f(X_y)

    out_samples[[length(out_samples)+1]] <- tibble(Y = Y, M = M, SampleID = r,
                                                   trt = data[[trt_name]])

    if(r %% print_interval == 0) cat(paste("Finishing iteration", r, "\n"))
  }

  out <- do.call(rbind, out_samples)
  return(out)

}
