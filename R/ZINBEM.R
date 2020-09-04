#' This function estimates the MLE of the parameters through Expected-Maximisation algorithm.
#'
#' @param Count Vector of count expression data of a gene over the cells.
#' @param group Vector which specifies the membership of the cells, i.e. two groups to be compared, corresponding to the cells in the count.
#' @param CellCluster Vector which specifies the cluster memberships of the cells, i.e. each entry represents memberships of the entries of the count.
#' @param CellAuxil Vector of cell level auxiliary information, corresponding to the entries in the counts, default is NULL.
#' @param maxit Maximum number of iterations for Expected-Maximization (EM) algorithm.
#' @param eps Convergennce criteria for EM algorithm.
#' @param muoffset Offset parameter for mean (mu) parameter, default is NULL.
#' @param phioffset Offset parameter for zero inflation (phi) parameter, default is NULL.
#' @param weights Observation wise weights for the cells, default is unity vector.
#'
#' @return A list containing the estimates of the parameters of a gene.
#'
#' @author Samarendra Das
#'
#' @import stats
#' @importFrom MASS glm.nb
#'
#' @export

ZINBEM <- function(Count, group, CellCluster, CellAuxil, maxit, eps, muoffset, phioffset, weights){

  loglikfun <- function(parms) {
    mu <- as.vector(exp(X %*% parms[1:kx] + offsetx))
    phi <- as.vector(plogis(Z %*% parms[(kx + 1):(kx + kz)] + offsetz))
    theta <- exp(parms[(kx + kz) + 1])
    loglik0 <- log(phi + exp(log(1 - phi) + suppressWarnings(dnbinom(0,
                                                                     size = theta, mu = mu, log = TRUE))))
    loglik1 <- log(1 - phi) + suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))
    loglik <- sum(weights[Y0] * loglik0[Y0]) + sum(weights[Y1] * loglik1[Y1])
    return(loglik)
  }

  ######design matrices
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
  ######iteration default parameters
  if(is.null(maxit))
    maxit <- 100
  if(is.null(eps))
    eps <- 1e-4
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
  #########starting values #family = MASS::negative.binomial(1)
  #model_count <- glm.fit(X, Y, family = poisson(link="log"), weights = weights, offset = offsetx)
  EM <- function(Y, Y0, Y1, X, Z, weights, offsetx, offsetz, maxit, eps){
    #cat ("\n" , paste("Finding starting values for EM algorithm..."))
    model_count <- suppressWarnings(glm.fit(X, Y, family = poisson(), weights = weights,
                                            offset = offsetx))
    model_zero <- suppressWarnings(glm.fit(Z, as.integer(Y0), weights = weights,
                                           family = binomial(link = "logit"), offset = offsetz))
    start <- list(count = model_count$coefficients, zero = model_zero$coefficients)
    start$theta <- 1
    start$count[is.na(start$count)]=0
    start$zero[is.na(start$zero)]=0
    mui <- model_count$fitted
    probi <- model_zero$fitted
    probi <- probi/(probi + (1 - probi) * dnbinom(0, size = start$theta, mu = mui))
    probi[Y1] <- 0
    ll_new <- loglikfun(c(start$count, start$zero, log(start$theta)))
    ll_old <- 2 * ll_new
    offset <- offsetx
    iter <- as.numeric()
    iter <- 0

    while (abs(ll_old - ll_new) > eps) {
      ll_old <- ll_new
      model_count <- suppressWarnings(glm.nb(Y ~ 0 +
                                               X + offset(offset), weights = weights*(1 -
                                                                                        probi), start = start$count, init.theta = start$theta))
      model_zero <- suppressWarnings(glm.fit(Z, probi,
                                             weights = weights, offset = offsetz, family = binomial(link = "logit"),
                                             start = start$zero))
      start <- list(count = model_count$coefficients,
                    zero = model_zero$coefficients, theta = model_count$theta)
      #print(c(start, it))
      mui <- model_count$fitted
      probi <- model_zero$fitted
      probi <- probi/(probi + (1 - probi) * dnbinom(0,
                                                    size = start$theta, mu = mui))
      probi[Y1] <- 0
      ll_new <- loglikfun(c(start$count, start$zero, log(start$theta)))
      iter <- iter + 1
      Convergence <- TRUE
      if (iter > maxit) {Convergence=FALSE; stop("Convergence not achieved in maxit")}
    }
    start$iter <- iter
    return(start)
  }

  start <- try( EM(Y = Y, Y0 = Y0, Y1= Y1, X = X, Z = Z, weights=weights, offsetx=offsetx, offsetz = offsetz, maxit =maxit, eps =eps), silent = TRUE)
  options(show.error.messages = TRUE)
  if('try-error' %in% class(start)){
    #cat("\n", paste("EM failed & Initiating optimization step..."))
    start <- ZINBoptim (Count=Y, group, CellCluster, CellAuxil, muoffset = offsetx,
                         phioffset=offsetz, weights)
    start$iter <- 0
  }

  #########count p-values
  #count_pval <- coef(summary(model_count))[,'Pr(>|z|)']
  #zero_pval <- coef(summary(model_zero))[,'Pr(>|z|)']
  ###
  mu <- exp(X %*%start$count  + offsetx)[, 1]
  mean.finite <- function(x) mean(x[is.finite(x)])
  pop_mean <- mean.finite(mu)
  phi <- plogis(Z %*% start$zero + offsetz)[, 1]
  pop_phi <- mean.finite(phi)
  theta <- start$theta
  ##########statistical testing for DE
  start$count0 <- start$count
  start$count0[2] <- 0
  mu0 <- exp(X %*%start$count0  + offsetx)[, 1]
  likH0 <- sum(dzinb(Y, size=start$theta, mu=mu0, rho = phi, log = TRUE))
  likH <- sum(dzinb(Y, size=start$theta, mu=mu, rho = phi, log = TRUE))
  lrt.stat <- -2 * (likH0 - likH)
  lrt.stat <- ifelse(lrt.stat<0, 0, lrt.stat)
  pval <- exp(pchisq(lrt.stat, df = 1, lower.tail = FALSE, log.p = TRUE))

  ############Differential Zero inflation
  start$zero0 <- start$zero
  #####under Null hypothesis
  start$zero0[2] <- 0
  phi0 <- plogis(Z %*%start$zero0  + offsetz)[, 1]
  phi.likH0 <- sum(dzinb(Y, size=start$theta, mu=mu, rho = phi0, log = TRUE))
  phi.likH <- sum(dzinb(Y, size=start$theta, mu=mu, rho = phi, log = TRUE))
  lrt.stat.phi <- -2 * (phi.likH0 - phi.likH)
  lrt.stat.phi <- ifelse(lrt.stat.phi < 0, 0, lrt.stat.phi)
  pval.phi <- exp(pchisq(lrt.stat.phi, df = 1, lower.tail = FALSE, log.p = TRUE))
  #######
  DE.stat <- c(lrt.stat, pval, lrt.stat.phi, pval.phi); names(DE.stat) <- c("Stat.DE","Pval.DE",
                                                                            "Stat.DZI", "Pval.DZI")
  #######Distribution characterization
  params <- c(pop_mean, pop_phi, theta)
  names(params) <- c("Mean", "ZeroInflation", "Dispersion")
  CountGroup <- start$count[1:2]; names(CountGroup) <- c("Intercept", "Group2")
  #CountGroup_p <- count_pval[1:2]; names(CountGroup_p) <- c("Intercept.p", "Group2.p")
  nclust <- max(as.numeric(levels(CellCluster)))
  countClust <- start$count[3: (nclust+1)]; names(countClust) <- paste("Cellcluster", 2:nclust)
  #countClust_p <- count_pval[3: (nclust+1)]; names(countClust_p) <- paste("pval.cluster", 2:nclust)
  ZeroGroup <- start$zero[1:2]; names(ZeroGroup) <- c("Intercept.0", "Group2.0")
  zeroClust <- start$zero[3: (nclust+1)]; names(zeroClust) <- paste("Cellcluster.0", 2:nclust)
  iter <- start$iter; names(iter) <- "#iteration"
  if(!is.null(CellAuxil)){
    nCellAux <- max(as.numeric(levels(cellTyp)))
    countCellAux <- start$count[(nclust+2):kx]; names(countCellAux) <- paste("CellAuxi", 2:nCellAux)
    #countCellAux_p <- count_pval[(nclust+2):kx]; names(countCellAux_p) <- paste("CellAuxi", 2:nCellAux)
    zeroCellAux <- start$zero[(nclust+2):kz]; names(zeroCellAux) <- paste("CellAuxi.0", 2:nCellAux)
    res <- c(params, CountGroup, countClust, countCellAux, ZeroGroup = ZeroGroup, zeroClust, zeroCellAux, iter)
    out <- list (res = res, DE.stat = DE.stat)
    remove(countCellAux, zeroCellAux, res, DE.stat, params, CountGroup, countClust,  ZeroGroup, zeroClust)
  } else {
    res <-  c(params, CountGroup, countClust, ZeroGroup, zeroClust, iter)
    out <- list( res=res, DE.stat = DE.stat)
    remove(res, DE.stat, params, CountGroup, countClust,  ZeroGroup, zeroClust,  iter)
  }
  return(out)
}
