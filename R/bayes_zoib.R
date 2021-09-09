# Save this file as `R/lm_stan.R`

#' Bayesian linear regression with Stan
#'
#' @export
#' @param x Numeric vector of input values.
#' @param y Numberic vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
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
