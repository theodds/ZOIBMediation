# ZOIBMediation

Mediation analysis for semicontinuous bounded data as performed in the paper
by Rene et al. (2021) Causal Mediation and Sensitivity Analysis for Mixed-Scale
Data. Additional details on the methodology can be found in that paper. The
package can be installed by running `devtools::install_github("theodds/ZOIBMediation")`.

We briefly illustrate some features of the package by reproducing the analysis
in that paper. First, we load the package and clean the JOBS II dataset
```{r}
library(ZOIBMediation)
library(dplyr)

jobs <- mediation::jobs |> dplyr::select(treat, econ_hard, depress1, sex,
                                  age, depress2, job_seek) |>
  dplyr::mutate(depress2 = (depress2 - 1) / 4, job_seek = (job_seek - 1) / 4,
         econ_hard = scale(econ_hard), depress1 = scale(depress1),
         age = scale(age))
```

Next, we fit the ZOIB model and perform g-computation to compute the mediation
effects under sequential ignorability. This might take a little bit of time.

```{r}
formula_m <- job_seek ~ treat + econ_hard + depress1 + sex + age
formula_y <- depress2 ~ job_seek + treat + econ_hard + depress1 + sex + age

jobs_fit <- bayes_zoib(
  formula_m = formula_m,
  formula_y = formula_y,
  df = jobs,
  cores = 4
)

mediated_jobs <- zoib_mediate(
  zoib_fit = jobs_fit, data = jobs, formula_m = job_seek ~ . - depress2,
  formula_y = depress2 ~ ., g_comp_thin = 1, print_interval = 10,
  med_name = "job_seek", trt_name = "treat", out_name = "depress2"
)
```

We can check mixing with the `traceplot_zoib` function and compute summary
statistics for the regression coefficients with `summary_zoib`.
```{r}
traceplot_zoib(fit = jobs_fit, formula = formula_m, data = jobs)
summary_zoib(jobs_fit, formula_y, formula_m, jobs)
```

We perform the sensitivity analysis to unmeasured confounding used the
`sensitivity_zoib` and `plot_sensitivity_zoib` functions:
```{r}
lambda_grid <- seq(from = -2, to = 2, by = 0.1)

sensitivity_zoib <- zoib_sensitivity(
  data = jobs, samples = jobs_fit, formula_m = formula_m, formula_y = formula_y,
  g_comp_thin = 4, print_interval = 1, med_name = "job_seek",
  trt_name = "treat", out_name = "depress2", lambda_grid = lambda_grid
)
plot_zoib_sensitivity(sensitivity_zoib)
```

Finally, `ppc_zoib` and `plot_ppc` perform posterior predictive checks.

```{r}
ppc_zoib <- zoib_ppc(
  samples = jobs_fit,
  thin = 100,
  data = jobs,
  formula_m = formula_m,
  formula_y = formula_y,
  med_name = "job_seek",
  trt_name = "treat"
)

plot_ppc(ppc_zoib, jobs, "treat", "job_seek", "depress2")
```


