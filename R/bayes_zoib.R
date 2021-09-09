#' Fit the zero-one inflated beta model for both the mediator and outcome
#'
#' For both the mediator and outcome we use (i) a logistic regression for the
#' probability of a 0, (ii) a logistic regression for the probability of 1
#' given non-zero, and (iii) a beta regression with logistic mean-link and
#' exponential precision link. That is:
#' \itemize{
#'   \item The event (M == 0) is modeled with a logistic regression
#'   \item The event (M == 1) is modeled with a conditiona logistic regression, where
#'   we work conditional on (M != 0)
#'   \item Given 0 < M < 1, we have M ~ Beta(mu * phi, (1 - mu) * phi) where
#'   mu = expit(X * beta_mu) and phi = exp(X * beta_phi)
#' }
#'
#' @export
#' @param formula_m Formula for the regression of the mediator on covariates and treatment .
#' @param formula_y Formula for the regression of the outcome on covariates, treatment, and the mediator
#' @param df A data.frame containing the covariates, mediator, outcome, and treatment
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling` which contains samples from the fitted ZOIB model.
#'
bayes_zoib <- function(formula_m, formula_y, df, ...) {

  frame_m <- model.frame(formula_m, data = df)
  frame_y <- model.frame(formula_y, data = df)

  X_m <- model.matrix(formula_m, data = df)
  X_y <- model.matrix(formula_y, data = df)
  y   <- model.response(frame_y)
  m   <- model.response(frame_m)

  zoib_data <- list(
    n = nrow(X_m),
    np_y = ncol(X_y),
    np_m = ncol(X_m),
    y = y,
    m = m,
    X_y = X_y,
    X_m = X_m
  )

  out <- rstan::sampling(stanmodels$bayes_zoib, data = zoib_data, ...)
  return(out)
}
