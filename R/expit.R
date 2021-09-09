#' Expit
#'
#' Returns the expit of x, i.e., (1 + exp(-x))^-1
#'
#' @param x a numeric vector
#'
#' @return Numeric vector exp(x) / (1 + exp(x))
expit <- function(x) {

  return(plogis(x))
}
