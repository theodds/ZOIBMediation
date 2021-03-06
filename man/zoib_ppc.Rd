% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zoib_pps.R
\name{zoib_ppc}
\alias{zoib_ppc}
\title{Posterior predictions from ZOIB marginals}
\usage{
zoib_ppc(
  samples,
  thin,
  data,
  formula_m,
  formula_y,
  med_name,
  trt_name,
  print_interval = 100
)
}
\arguments{
\item{samples}{Model fit with the bayes_zoib function}

\item{thin}{Do we want to sample on all iterations (thin = 1) or half (thin = 2), etc.}

\item{data}{The data used to fit the model}

\item{formula_m}{Formula used to fit the mediator model}

\item{formula_y}{Formula used to fit the outcome model}

\item{med_name}{Name of the mediator}

\item{trt_name}{Name of the treatment}

\item{print_interval}{Number of iterations before printing an update on progress}
}
\value{
Samples from the posterior predictive distribution, stored as a data.frame
giving the subject id, iteration, outcome, and mediator.
}
\description{
Draws from the posterior predictive distribution of the mediator and outcome
}
