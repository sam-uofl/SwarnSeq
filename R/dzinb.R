#' # This function computes the distribution of the zero inflated negative binomial distribution.
#'
#' @param x Vector of (non-negative integer) quantiles.
#' @param mu mean parameter
#' @param size size parameter.
#' @param rho Zero inflation parameter.
#' @param log Logical variable either TRUE/FALSE; if TRUE, probabilities p are given as log(p).
#'
#' @return Returns the distribution function for the zero inflated negative binomial distribution with mean, size, zero inflation parameters.
#'
#' @author Samarendra Das
#'
#' @importFrom stats dnbinom
#'
#' @examples
#'
#' x = c(0, 1, 4)
#' mu <- 0.5; size <- 0.23; rho <- 0.32
#' p = dzinb(x, mu, size, rho, log=TRUE)
#'
#' @export
#'
dzinb <- function (x, size, mu, rho, log)
{
  if(sum(is.na(x)) > 0)
    stop("NA detected in x");gc();
  if(sum(x < 0) > 0)
    stop("Negative values detected in x");
  gc()
  if (size < 0 || mu < 0 || rho < 0){
    warning ("Parameters, i.e. mu, size, phi, must be positive")
    return(NaN); gc()
  }
  if (rho < 0 || rho > 1)
    warning ("Probability of structural zero's must be in [0, 1]"); gc()
  if (!is.logical(log) || length(log) != 1)
    stop("Bad input for argument 'log', must be either TRUE or FALSE"); gc()
  out <- rho * (x==0) + (1 - rho) * dnbinom (x, size = size, mu = mu, log = FALSE)
  if(missing(log))
    log <- TRUE
  if (log==TRUE)
    log(out)
  else if (log==FALSE)
    return(out)
}

