#' Traceplots for ZOIB fit
#'
#' Returns a traceplot for the coefficients of a fitted ZOIB model
#'
#' @param fit Model fit with bayes_zoib
#' @param formula Formula used to fit the model (either for the mediator or outcome, depending on the response type)
#' @param data Data used to fit the model
#' @param param The model component we want to look at; can be "alpha", "gamma",
#' "mu", or "phi" where alpha corresponds to the mass at 0, gamma to the mass at
#'  1, and (mu,phi) to the mean/precision of the beta model.
#' @param response "mediator" or "outcome"?
#'
#' @return A ggplot object showing traceplots.
#' @export
traceplot_zoib <- function(fit, formula, data,
                           param = "alpha", response = "mediator") {

  Val <- iteration <- chain <- NULL

  vnames <- colnames(model.matrix(formula, data = data))
  stopifnot(response %in% c("mediator", "outcome"))

  num_chains <- length(fit@stan_args)
  num_after_burn <- fit@stan_args[[1]][["iter"]] - fit@stan_args[[1]][["warmup"]]
  num_save <- num_after_burn / fit@stan_args[[1]][["thin"]]

  stopifnot(param %in% c("alpha", "gamma", "mu", "phi"))
  col_idx <- 1
  col_idx <- ifelse(param == "gamma", 2, col_idx)
  col_idx <- ifelse(param == "mu", 3, col_idx)
  col_idx <- ifelse(param == "phi", 4, col_idx)

  sample_mat <- as.matrix(fit)
  sample_df_med_alpha  <- sample_mat |> as.data.frame() |>
    select(starts_with("beta")) |>
    select(contains(response)) |>
    select(contains(paste0(col_idx,"]"))) |>
    mutate(iteration = rep(1:num_save, num_chains), chain = rep(1:num_chains, each = num_save)) %>%
    pivot_longer(cols = starts_with("beta"), names_to = "Param", values_to = "Val")

  new_params <- sample_df_med_alpha$Param
  old_names <- unique(sample_df_med_alpha$Param)
  for(i in 1:length(vnames)) {
    new_params[new_params == old_names[i]] <- vnames[i]
  }
  sample_df_med_alpha_2 <- sample_df_med_alpha |> mutate(Param = new_params)

  my_plot <- ggplot(sample_df_med_alpha_2,
                    aes(y = Val, x = iteration, color = factor(chain))) +
    geom_line(alpha = .5) +
    facet_wrap(~Param, ncol = 5, scale = "free_y") +
    # scale_color_viridis_d() +
    theme_bw() + theme(legend.position = "none") + xlab("Iteration") + ylab("")

  return(my_plot)

}
