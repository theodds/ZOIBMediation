#' logit
#'
#' Returns the logit of x
#'
#' @param x numeric vector between 0 and 1
#'
#' @return Numeric vector log(x) - log(1 - x)
#' @export
logit <- function(x) {
  qlogis(x)
}
