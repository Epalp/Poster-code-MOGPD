#'  Marshall Olkins Generalised Pareto Distribution family
#'
#' @description
#' The function \code{MOGPD()} defines the Marshall Olkins Generalised Pareto Distribution, a three parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "identity" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#' @param nu.link defines the nu.link, with "log" link as the default for the nu.
#'
#' @details
#' The Generalised exponential-Gaussian with parameters \code{mu}, \code{sigma} and \code{nu}
#' has density given by
#'
#'
#'\eqn{f(y | \mu, \sigma, \nu) = \frac{\frac{\mu}{\nu}*(1+\frac{\sigma*y}{\nu})^{(\frac{-1}{\sigma}-1)}}{(1-(1-\mu)*(1+\frac{\sigma*y}{\nu})^{(\frac{-1}{\sigma})})^2}}
#'
#' @export
MOGPD <- function (mu.link="log", sigma.link="log", nu.link="log") {
  mstats <- checklink("mu.link", "Marshall Olkins Generalised Pareto Distribution ",
                      substitute(mu.link), c("log", "own"))
  dstats <- checklink("sigma.link", "Marshall Olkins Generalised Pareto Distribution",
                      substitute(sigma.link), c("log", "own"))
  vstats <- checklink("nu.link", "Marshall Olkins Generalised Pareto Distribution",
                      substitute(nu.link), c("log", "own"))

  structure(list(family=c("MOGPD", "Marshall Olkins Generalised Pareto Distribution"),
                 parameters=list(mu=TRUE, sigma=TRUE, nu=TRUE),
                 nopar=3,
                 type="Continuous",

                 mu.link = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),
                 nu.link = as.character(substitute(nu.link)),

                 mu.linkfun = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,
                 nu.linkfun = vstats$linkfun,

                 mu.linkinv = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,
                 nu.linkinv = vstats$linkinv,

                 mu.dr = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,
                 nu.dr = vstats$mu.eta,

                 # Primeras derivadas ---------------------------------
                 dldm = function(y, mu, sigma, nu) {
                   p1<-mu+(nu*y/sigma+1)^(1/sigma)-1
                   p2<-mu*(mu+(nu*y/sigma+1)^(1/sigma)-1)
                   dldm<-p1/p2
                   dldm
                 },

                 dldd = function(y, mu, sigma, nu) {
                   p1<-log(nu*y/sigma+1)
                   p2<-sigma^2
                   dldd<-p1/p2
                   dldd
                 },

                 dldv = function(y, mu, sigma, nu) {
                   p1<-mu*(sigma-2)+sigma*((nu*y/sigma+1)^(1/sigma)-1)+2
                   p2<-sigma*(mu+(nu*y/sigma+1)^(1/sigma)-1)
                   p3<-sigma+nu*y
                   dldv<-(-1/sigma-(p1/p2))/p3
                   dldv
                 },

                 # Segundas derivadas ---------------------------------
                 d2ldm2 = function(y, mu, sigma, nu) {
                   p1<-mu^2-2*mu*((sigma*y/nu+1)^(1/sigma)-1)-((sigma*y/nu+1)^(1/sigma)-1)^2
                   p2<-mu^2*(mu+(sigma*y/nu+1)^(1/sigma)-1)^2
                   d2ldm2<-p1/p2
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma, nu) {

                   A<-mu+(nu*y/sigma+1)^(1/sigma)-1
                   B<-mu*(mu+(nu*y/sigma+1)^(1/sigma)-1)
                   dldm<-A/B

                   C<-log(nu*y/sigma+1)
                   D<-sigma^2
                   dldd<-C/D

                   d2ldmdd <- -dldm * dldd
                   d2ldmdd
                 },

                 d2ldmdv = function(y, mu, sigma, nu) {

                   A<-mu+(nu*y/sigma+1)^(1/sigma)-1
                   B<-mu*(mu+(nu*y/sigma+1)^(1/sigma)-1)
                   dldm<-A/B

                   C<-mu*(sigma-2)+sigma*((nu*y/sigma+1)^(1/sigma)-1)+2
                   D<-sigma*(mu+(nu*y/sigma+1)^(1/sigma)-1)
                   E<-sigma+nu*y
                   dldv<-(-1/sigma-(C/D))/E

                   d2ldmdd <- -dldm * dldv
                   d2ldmdd
                 },

                 d2ldd2 = function(y, mu, sigma, nu) {
                   p1<--2*sigma*(sigma*y^2*(sigma*(-mu^2+2*mu+(sigma*y/nu+1)^(2/sigma)-1))-2*(mu-1)*(sigma*y/nu+1)^(1/sigma))
                   p2<-nu^2*(-mu^2+2*mu+(sigma*y/nu+1)^(2/sigma)-1)
                   p3<-2*nu*y*(sigma*(-mu^2+2*mu+(sigma*y/nu+1)^(2/sigma)-1)-(mu-1)*(sigma*y/nu)^(1/sigma))
                   p4<-log(sigma*y/nu+1)+sigma^2*y*(y*(3*sigma*(-mu^2+2*mu+(sigma*y/nu+1)^(2/sigma)-1)+sigma^2*(mu+(sigma*y/nu+1)^(1/sigma)-1)^2-2*(mu-1)*(sigma*y/nu+1)^(1/sigma))+2*nu*(-mu^2+2*mu+(sigma*y/nu+1)^(2/sigma)-1))
                   p5<--2*(mu-1)*(nu+sigma*y)^2*(sigma*y/nu+1)^(1/sigma)*log(sigma*y/nu+1)
                   p6<-(sigma^4*(nu+sigma*y)^2*((sigma*y/nu+1)^(1/sigma)-1)^2)
                   d2ldd2<-(p1+p2+p3+p4+p5)/p6
                   d2ldd2
                 },

                 d2ldddv = function(y, mu, sigma, nu) {

                   A<-log(nu*y/sigma+1)
                   B<-sigma^2
                   dldd<-A/B

                   C<-mu*(sigma-2)+sigma*((nu*y/sigma+1)^(1/sigma)-1)+2
                   D<-sigma*(mu+(nu*y/sigma+1)^(1/sigma)-1)
                   E<-sigma+nu*y
                   dldv<-(-1/sigma-(C/D))/E

                   d2ldddv <- -dldd * dldv
                   d2ldddv
                 },

                 d2ldv2 = function(y, mu, sigma, nu) {
                   p1<-y^2*(sigma*(mu^2-2*mu-(sigma*y/nu+1)^(2/sigma)+1)-2*(mu-1)*(sigma*y/nu+1)^(1/sigma))
                   p2<-nu^2*(mu+(sigma*y/nu+1)^(1/sigma)-1)^2-2*nu*y*(-mu^2+2*mu+(sigma*y/nu+1)^(2/sigma)-1)
                   p3<-(nu^2*(nu+sigma*y)^2*(mu+(sigma*y/nu+1)^(1/sigma)-1)^2)
                   d2ldv2<-(p1+p2)/p3
                   d2ldv2
                 },

                 G.dev.incr = function(y, mu, sigma, nu, ...) -2*dMOGPD(y, mu, sigma, nu, log=TRUE),
                 rqres = expression(rqres(pfun="pMOGPD", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),

                 mu.initial = expression(mu       <- rep(9, length(y))),
                 sigma.initial = expression(sigma <- rep(1, length(y))),
                 nu.initial = expression(nu       <- rep(1, length(y))),

                 mu.valid = function(mu)       all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),
                 nu.valid = function(nu)       all(nu > 0),

                 y.valid = function(y) all(y > 0)
  ),
  class=c("gamlss.family", "family"))
}
