sim_zoib_copula <- function(u, alpha, gamma, mu, phi) {
  u_prime <- (u - alpha) / ((1 - alpha) * (1 - gamma))
  N <- length(u)
  out <- numeric(N)
  out[u < alpha] <- 0
  out[u > 1 - (1 - alpha) * gamma] <- 1
  out <- ifelse(0 < u_prime & u_prime < 1,
                qbeta(u_prime, mu * phi, (1 - mu) * phi), out)
  return(out)

}

sim_zoib_copula_df <- function(u, X, beta_alpha, beta_gamma, beta_mu, beta_phi) {
  alpha <- expit(X %*% beta_alpha) |> as.numeric()
  gamma <- expit(X %*% beta_gamma) |> as.numeric()
  mu    <- expit(X %*% beta_mu)    |> as.numeric()
  phi   <- exp(  X %*% beta_phi)   |> as.numeric()
  sim_zoib_copula(u = u, alpha = alpha, gamma = gamma, mu = mu, phi = phi)
}

mean_zoib <- function(alpha, gamma, mu) {
  (1 - alpha) * gamma + (1 - alpha) * (1 - gamma) * mu
}

mean_zoib_df <- function(X, beta_alpha, beta_gamma, beta_mu, beta_phi) {
  alpha <- expit(X %*% beta_alpha) |> as.numeric()
  gamma <- expit(X %*% beta_gamma) |> as.numeric()
  mu    <- expit(X %*% beta_mu)    |> as.numeric()
  phi   <- exp(  X %*% beta_phi)   |> as.numeric()
  mean_zoib(alpha = alpha, gamma = gamma, mu = mu)
}
