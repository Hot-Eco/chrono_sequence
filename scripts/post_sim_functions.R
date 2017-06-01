#' Draws random values from a mulitvariate normal distribution
#' 
#' @param n Number of values required.
#' @param mu Vector of means.
#' @param vc Variance-covariance matrix
#' 
#' @return A matrix of random values with n rows and p columns where
#'   \code{p = length(mu)}.
#' 
rmvn <- function(n, mu, vc) {
  L <- mgcv::mroot(vc)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m * n), m, n))
}


#' Generates random parameter values from a fitted model
#' 
#' @param n Number of values required.
#' @param model The fitted model (e.g. a glm or gam object).
#' 
#' @return A matrix of random values with n rows and p columns where
#'   \code{p = length(coef(model))}.
#' 
rparams <- function(n, model) {
  rmvn(n, coef(model), vcov(model))
}


#' Generates posterior predictions from a fitted model
#' 
#' @param n Number of predictions required.
#' 
#' @param model A fitted model (e.g. a glm or gam object).
#' 
#' @param newdata A data frame defining the grid over which
#'   to generate predictions (e.g. as would be used with the
#'   standard \code{predict} function).
#'   
#' @param prefix A character string to use as a prefix for
#'   column labels identifying simulation replicates.
#'   
#' @param join If \code{TRUE}, append the simulated data as
#'   extra columns to \code{newdata}; otherwise return just
#'   the simulated data alone.
#' 
#' @return A data.frame of simulated data with, if \code{join}
#'   was \code{TRUE}, the columns from \code{newdata} pre-pended.
#' 
postSim <- function(n, model, newdata, 
                    prefix = "sim.",
                    join = TRUE) {
  
  lp <- predict(model, newdata = newdata, type = "lpmatrix")
  rb <- rparams(n, model)
  sim <- lp %*% t(rb)
  
  prefix <- prefix[1] # just in case
  colnames(sim) <- paste(prefix, 1:ncol(sim), sep="")
  
  if (join) cbind(dat.predict, sim)
  else sim
}

