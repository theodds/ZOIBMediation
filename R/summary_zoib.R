summary_zoib <- function(fit, formula_y, formula_m, data) {

  vnames_y <- paste0(colnames(model.matrix(formula_y, data = data)))
  vnames_m <- paste0(colnames(model.matrix(formula_m, data = data)))

  posterior_summary <- fit |>
    rstan::summary(pars = c("beta_mediator", "beta_outcome"), probs = c(0.025, 0.975)) |>
    pluck("summary") |>
    as.data.frame() |>
    mutate(Model = c(rep(c("alpha", "gamma", "mu", "phi"), each = length(vnames_m)),
                     rep(c("alpha", "gamma", "mu", "phi"), each = length(vnames_y))),
           Response = c(rep("Mediator", 4 * length(vnames_m)), rep("Outcome", 4 * length(vnames_y))),
           Variable = c(rep(vnames_m, 4), rep(vnames_y,4))) |>
    relocate("Variable", "Model", "Response")
  rownames(posterior_summary) <- NULL

  return(posterior_summary)
}
