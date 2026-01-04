#' Modified from EpiSki (https://github.com/alicia-gill/EpiSky/tree/master)
#' 
#' Skellam Log-Likelihood
#'
#' Approximates the probability of a change in prevalence given the birth and death rate.
#'
#' @param new_x new prevalence.
#' @param old_x old prevalence.
#' @param birth_rate birth rate of the epidemic.
#' @param death_rate death rate of the epidemic.
#' @param log logical; if TRUE, log-likelihood is given.
#'
#' @return log-likelihood.
#' @export
#'
dskellam_sir2 <- function(new_x, old_x, birth_rate, death_rate, log=T) {
  
  new_x <- as.numeric(new_x)
  old_x <- as.numeric(old_x)
  nu <- abs(new_x - old_x)
  x <- 2 * old_x * sqrt(birth_rate * death_rate)
  length <- length(nu)
  logI <- rep(NA, length)
  nu_cutoff <- 10
  x_cutoff <- 10000
  set0 <- (1:length)[nu <= nu_cutoff & x <= x_cutoff] 
  set1 <- (1:length)[nu <= nu_cutoff & x > x_cutoff]
  set2 <- (1:length)[nu > nu_cutoff]
  
  if (length(birth_rate) > 1) {
    logI[set0] <- log(besselI(x = 2 * old_x[set0] * sqrt(birth_rate[set0] * death_rate), nu = nu[set0], expon.scaled = T))
    logI[set1] <- Bessel::besselIasym(x = 2 * old_x[set1] * sqrt(birth_rate[set1] * death_rate), nu = nu[set1], k.max = 20, expon.scaled = T, log = T)
    logI[set2] <- Bessel::besselI.nuAsym(x = 2 * old_x[set2] * sqrt(birth_rate[set2] * death_rate), nu = nu[set2], k.max = 5, expon.scaled = T, log = T)
  } else {
    logI[set0] <- log(besselI(x = 2 * old_x[set0] * sqrt(birth_rate * death_rate), nu = nu[set0], expon.scaled = T))
    logI[set1] <- Bessel::besselIasym(x = 2 * old_x[set1] * sqrt(birth_rate * death_rate), nu = nu[set1], k.max = 20, expon.scaled = T, log = T)
    logI[set2] <- Bessel::besselI.nuAsym(x = 2 * old_x[set2] * sqrt(birth_rate * death_rate), nu = nu[set2], k.max = 5, expon.scaled = T, log = T)
  }
  
  logl <- ((new_x - old_x)/2 * (log(birth_rate) - log(death_rate))) -
    (old_x * (birth_rate + death_rate)) +
    logI +
    (2 * old_x * sqrt(birth_rate * death_rate))
  
  if (log==T) {
    return(logl)
  } else {
    return(exp(logl))
  }
}