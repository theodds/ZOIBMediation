#' Posterior predictive checks of marginals
#'
#' Makes plots assessing the fit of the posterior predictive distribution to the
#' data. Returns two figures as ggplot objects: the first gives the posterior
#' predictive distribution of a kernel density estimate to the continuous part
#' of the outcome/mediator stratified by treatment, the second gives posterior
#' samples of the proportions of 0s and 1s for the outcome/mediator stratified
#' by treatment. For both figures, the statistics computed on the original data
#' are provided as a baseline.
#'
#' @param ppc An objected obtained from the fitted model using the zoib_ppc function
#' @param data The original dataset
#' @param trt_name Name of the treatment variable in the original dataset
#' @param med_name Name of the mediator variable in the original dataset
#' @param out_name Name of the outcome variable in the original dataset
#'
#' @return ggplots giving the posterior probability check for the continuous
#' and discrete parts of the data
#' @export
plot_ppc <- function(ppc, data, trt_name, med_name, out_name) {

  Variable <- Y <- M <- trt <- Value <- SampleID <- Proportion <- NULL

  df_plot_1 <- ppc |>
    pivot_longer(cols = c("Y", "M"), names_to = "Variable", values_to = "Value") |>
    mutate(Variable = ifelse(Variable == "Y", out_name, med_name),
           trt = ifelse(trt == 1, "Treated", "Untreated"))
  df_plot_2 <- data
  df_plot_2$Y <- data[[out_name]]
  df_plot_2$M <- data[[med_name]]
  df_plot_2$trt <- ifelse(data[[trt_name]] == 1, "Treated", "Untreated")
  df_plot_2 <- df_plot_2 |> select(Y, M, trt) |>
    pivot_longer(cols = c("Y", "M"), names_to = "Variable", values_to = "Value") |>
    mutate(Variable = ifelse(Variable == "Y", out_name, med_name))

  p_1 <- ggplot() + geom_density(
    data = df_plot_1 |> filter(Value != 0, Value != 1),
    mapping = aes(x = Value, group = SampleID),
    color = scales::muted("red", l = 70, c = 70), alpha = 0.4) +
    facet_grid(Variable ~ factor(trt)) + theme_bw()

  p_2 <- geom_density(data = df_plot_2 |> filter(Value != 0, Value != 1),
                      mapping = aes(x = Value))

  density_plot <- p_1 + p_2 + xlab("") + ylab("Density")

  df_plot_3 <- df_plot_1 |>  group_by(Variable, SampleID, trt) |>
    summarise(`0` = mean(Value == 0), `1` = mean(Value == 1)) |>
    pivot_longer(cols = c("0", "1"),
                 names_to = "Value",
                 values_to = "Proportion")

  df_plot_4 <- df_plot_2 |> group_by(Variable, trt) |>
    summarise(`0` = mean(Value == 0), `1` = mean(Value == 1)) |>
    pivot_longer(cols = c("0", "1"),
                 names_to = "Value",
                 values_to = "Proportion")
  mass_plot <- ggplot(df_plot_3, aes(x = Proportion, fill = Value)) +
    geom_histogram() +
    facet_grid(Variable ~ trt) +
    geom_vline(mapping = aes(xintercept = Proportion), data = df_plot_4,
               alpha = 1, size = 1, linetype = 2, color = "pink") +
    theme_bw() + xlab("Proportion of Observations") + ylab("Frequency") + scale_fill_viridis_d(option = "E")

  return(list(density_plot = density_plot, mass_plot = mass_plot))

}
