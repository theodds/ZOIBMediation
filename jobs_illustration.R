library(ZOIBMediation)
library(zeallot)

## Load data ----

jobs <- mediation::jobs |> dplyr::select(treat, econ_hard, depress1, sex,
                                  age, depress2, job_seek) |>
  dplyr::mutate(depress2 = (depress2 - 1) / 4, job_seek = (job_seek - 1) / 4,
         econ_hard = scale(econ_hard), depress1 = scale(depress1),
         age = scale(age))

## Fit model and do mediation ----

formula_m <- job_seek ~ treat + econ_hard + depress1 + sex + age
formula_y <- depress2 ~ job_seek + treat + econ_hard + depress1 + sex + age

# jobs_fit <- bayes_zoib(
#   formula_m = formula_m,
#   formula_y = formula_y,
#   df = jobs,
#   cores = 4
# )
#
# mediated_jobs <- zoib_mediate(
#   zoib_fit = jobs_fit, data = jobs, formula_m = job_seek ~ . - depress2,
#   formula_y = depress2 ~ ., g_comp_iters = 1:4000, print_interval = 10,
#   med_name = "job_seek", trt_name = "treat", out_name = "depress2"
# )

# saveRDS(list(jobs_fit = jobs_fit, mediated_jobs = mediated_jobs),
# file = "mediated_zoib.rds")

c(jobs_fit, mediated_jobs) %<-% readRDS("mediated_zoib.rds")

## Check mixing ----

traceplot_zoib(fit = jobs_fit, formula = formula_m, data = jobs)

## Examine coefficients ----

summary_zoib(jobs_fit, formula_y, formula_m, jobs)

## Summarize mediation effects ----


