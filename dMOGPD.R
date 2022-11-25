#' The Marshall Olkins Generalised Pareto Distribution
#'
#' Marshall Olkins Generalised Pareto Distribution family
#'
#' @description
#' The function \code{MOGPD()} defines the Marshall Olkins Generalised Pareto Distribution, a three parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param y	vector of quantiles.
#' @param n number of observations.
#' @param mu parameter.
#' @param sigma parameter.
#' @param nu parameter.
#' @param log if TRUE, probabilities p are given as log(p).
#'
#' @details
#' The  Marshall Olkins Generalised Pareto Distribution with parameters \code{mu}, \code{sigma} and \code{nu}
#' has density given by
#'
#'\eqn{f(y | \mu, \sigma, \nu) = \frac{\frac{\mu}{\nu}*(1+\frac{\sigma*y}{\nu})^{(\frac{-1}{\sigma}-1)}}{(1-(1-\mu)*(1+\frac{\sigma*y}{\nu})^{(\frac{-1}{\sigma})})^2}}
#'
#'
#' @export
dMOGPD <- function(y,mu,sigma,nu,log=FALSE){
  if (any(mu <= 0))
    stop("parameter mu has to be positive!")
  if (any(sigma <= 0))
    stop("parameter sigma has to be positive!")
  if (any(nu <= 0))
    stop("parameter nu has to be positive!")
  res <- ((mu/nu)*(1+((sigma*y)/nu))^(-1/sigma-1))/(1-(1-mu)*(1+((sigma*y)/nu))^(-1/sigma))^2

  if(log)
    res<-log(res)
  return(res)
}
#' @importFrom stats uniroot
#' @export
#' @rdname dMOGPD

pMOGPD<-function(q, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  cdf <- 1-(1+((sigma*q)/nu))^(-1/sigma)
cdf
}

qMOGPD <- function(p, mu, sigma, nu, lower.tail=TRUE, log.p=FALSE){
  if (any(mu<= 0))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  fda1 <- function(q, mu, sigma, nu, p) {
    pMOGPD(q, mu, sigma, nu) - p
  }

  r_de_la_funcion <- function(mu, sigma, nu, p) {
    uniroot(fda1, interval=c(0, 1e+30), mu, sigma, nu, p)$root
  }
  r_de_la_funcion <- Vectorize(r_de_la_funcion)
  q <- r_de_la_funcion(mu, sigma, nu, p)
  q
}
#' @importFrom stats runif
#' @importFrom stats uniroot
#' @export
#' @rdname dMOGPD
rMOGPD <- function(n, mu, sigma, nu){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))

  n <- ceiling(n)
  p <- runif(n)
  r <- qMOGPD(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dMOGPD
rMOGPD(10,9,0.25,0.25)
