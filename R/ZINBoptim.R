#' #This function optimises the likelihood function for a gene, supplies starting values for EM algorithm.
#'
#' @param Count Vector of count expression data of a gene over the cells.
#' @param group Vector which specifies the membership of the cells, i.e. two groups to be compared, corresponding to the cells in the count.
#' @param CellCluster Vector which specifies the cluster memberships of the cells, i.e. each entry represents memberships of the entries of the count.
#' @param CellAuxil Vector of cell level auxiliary information, corresponding to the entries in the counts, default is NULL.
#' @param muoffset Offset parameter for mean (mu) parameter, default is NULL.
#' @param phioffset Offset parameter for zero inflation (phi) parameter, default is NULL.
#' @param weights Observation wise weights for the cells, default is unity vector.
#'
#' @return A list containing the estimates of the parameters of a gene.
#'
#' @author Samarendra Das
#'
#' @import stats
#'
#' @export
#'
ZINBoptim <- function(Count, group, CellCluster, CellAuxil, muoffset, phioffset, weights){
  gradNegBin <- function(parms) {
    eta <- as.vector(X %*% parms[1:kx] + offsetx)
    mu <- exp(eta)
    etaz <- as.vector(Z %*% parms[(kx + 1):(kx + kz)] + offsetz)
    muz <- plogis(etaz)
    theta <- exp(parms[(kx + kz) + 1])
    clogdens0 <- dnbinom(0, size = theta, mu = mu, log = TRUE)
    dens0 <- muz * (1 - as.numeric(Y1)) + exp(log(1 - muz) +
                                                clogdens0)
    wres_count <- ifelse(Y1, Y - mu * (Y + theta)/(mu + theta),
                         -exp(-log(dens0) + log(1 - muz) + clogdens0 + log(theta) -
                                log(mu + theta) + log(mu)))
    wres_zero <- ifelse(Y1, -1/(1 - muz) * make.link("logit")$mu.eta(etaz),
                        (make.link("logit")$mu.eta(etaz) - exp(clogdens0) * make.link("logit")$mu.eta(etaz))/dens0)
    wres_theta <- theta * ifelse(Y1, digamma(Y + theta) -
                                   digamma(theta) + log(theta) - log(mu + theta) + 1 -
                                   (Y + theta)/(mu + theta), exp(-log(dens0) + log(1 -
                                                                                     muz) + clogdens0) * (log(theta) - log(mu + theta) +
                                                                                                            1 - theta/(mu + theta)))
    colSums(cbind(wres_count * weights * X, wres_zero * weights *
                    Z, wres_theta))
  }
  #####likelihood
  loglikfun <- function(parms) {
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(plogis(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
    theta <- exp(parms[(kx + kz) + 1])
    loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0,
                                                                     size = theta, mu = mu, log = TRUE))))
    loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
    loglik <- -sum(weights[Y0] * loglik0[Y0]) - sum(weights[Y1] * loglik1[Y1])
    return(loglik)
  }
  group <- as.factor(group)
  CellCluster <- as.factor(CellCluster)
  if(is.null(CellAuxil)) {
    X <- model.matrix(~ group + CellCluster)
    Z <- model.matrix(~ group + CellCluster)
  }else {
    cellTyp <- as.factor(CellAuxil)
    X <- model.matrix(~ group + CellCluster + cellTyp)
    Z <- model.matrix(~ group + CellCluster + cellTyp)
  }
  Y <- Count
  remove(Count)
  n <- length(Y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  Y0 <- Y <= 0
  Y1 <- Y > 0
  if (is.null(weights))
    weights <- 1
  if (length(weights) == 1)
    weights <- rep.int(weights, n)
  weights <- as.vector(weights)

  #######OffsetX
  if (is.null(muoffset))
    muoffset <- 0
  offsetx <- muoffset
  if (length(offsetx) == 1)
    offsetx <- rep.int(offsetx, n)
  offsetx <- as.vector(offsetx)
  #######Offsetz
  if (is.null(phioffset))
    phioffset <- 0
  offsetz <- phioffset
  if (length(offsetz) == 1)
    offsetz <- rep.int(offsetz, n)
  offsetz <- as.vector(offsetz)
  remove(muoffset, phioffset)
  #########starting values for optim function
  #message("Finding starting values for optim function ...")
  model_count <- suppressWarnings(glm.fit(X, as.numeric(Y), family = poisson(), weights = weights,
                                          offset = offsetx))
  model_zero <- suppressWarnings(glm.fit(Z, as.integer(Y0), weights = weights,
                                         family = binomial(link = "logit"), offset = offsetz))
  start <- list(count = model_count$coefficients, zero = model_zero$coefficients)
  start$theta <- 1
  start$count[is.na(start$count)]=0
  start$zero[is.na(start$zero)]=0
  fit <- optim(fn = loglikfun, gr = gradNegBin, par = c(start$count,
                                                        start$zero, log(start$theta)),
               method = "BFGS", hessian = FALSE, control = list())

  start <- list()
  start$count <- fit$par[1:kx]
  names(start$count) <- names(start$count) <- colnames(X)
  start$zero <- fit$par[(kx + 1):(kx + kz)]
  names(start$zero) <- names(start$zero) <- colnames(Z)
  np <- kx + kz + 1
  start$theta <- exp(fit$par[np])
  #names(start$theta) <- "theta"
  return(start)
}
