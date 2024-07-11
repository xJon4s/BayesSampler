#' Metropolis-Hastings Sampling
#'
#' @param target_density Function for the target density to sample from.
#' @param proposal_sd Standard deviation of the proposal distribution.
#' @param start Initial value for the chain.
#' @param n_samples Number of samples to generate.
#' @return A numeric vector of samples.
#' @examples
#' target <- function(x) { dnorm(x, mean = 0, sd = 1) }
#' samples <- metropolis_hastings(target, proposal_sd = 1, start = 0, n_samples = 1000)
#' @export
metropolis_hastings <- function(target_density, proposal_sd, start, n_samples) {
  samples <- numeric(n_samples)
  samples[1] <- start
  for (i in 2:n_samples) {
    proposal <- rnorm(1, mean = samples[i - 1], sd = proposal_sd)
    accept_prob <- min(1, target_density(proposal) / target_density(samples[i - 1]))
    if (runif(1) < accept_prob) {
      samples[i] <- proposal
    } else {
      samples[i] <- samples[i - 1]
    }
  }
  return(samples)
}

#' Gibbs Sampling for Bivariate Normal Distribution
#'
#' @param mu Vector of means for the bivariate normal distribution.
#' @param sigma Covariance matrix for the bivariate normal distribution.
#' @param start Initial values for the chain.
#' @param n_samples Number of samples to generate.
#' @return A matrix of samples.
#' @examples
#' mu <- c(0, 0)
#' sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' samples <- gibbs_sampler(mu, sigma, start = c(0, 0), n_samples = 1000)
#' @export
gibbs_sampler <- function(mu, sigma, start, n_samples) {
  samples <- matrix(NA, n_samples, 2)
  samples[1, ] <- start
  for (i in 2:n_samples) {
    samples[i, 1] <- rnorm(1, mean = mu[1] + sigma[1, 2] / sigma[2, 2] * (samples[i - 1, 2] - mu[2]), sd = sqrt(sigma[1, 1] - sigma[1, 2]^2 / sigma[2, 2]))
    samples[i, 2] <- rnorm(1, mean = mu[2] + sigma[2, 1] / sigma[1, 1] * (samples[i, 1] - mu[1]), sd = sqrt(sigma[2, 2] - sigma[2, 1]^2 / sigma[1, 1]))
  }
  return(samples)
}

