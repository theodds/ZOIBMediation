#' Plot sensitivity analysis for ZOIB model
#'
#' @param sensitivity Objected obtained from the zoib_sensitivity function
#'
#' @return Plot of the sensitivity analysis for different values of the sensitivity
#' parameter lambda, as defined by the zoib_sensitivity function.
#' @export
plot_zoib_sensitivity <- function(sensitivity) {

  ParamName <- lambda <- Param <- LCL <- UCL <- Est <- NULL

  sensitivity |> group_by(ParamName, lambda) |>
    filter(ParamName != "EM0", ParamName != "EM1", ParamName != "tau") |>
    summarise(Est = mean(Param),
              LCL = quantile(Param, 0.025),
              UCL = quantile(Param, 0.975)) |>
    ggplot(aes(x = lambda)) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL),
                fill = "darkseagreen4",
                alpha = 0.3) +
    geom_line(aes(y = Est), linetype = 2) +
    facet_wrap(~ ParamName, nrow = 2) + theme_bw() +
    xlab("lambda") + ylab("Estimate") ->
    sensitivity_plot

  return(sensitivity_plot)
}
