#' This function decides the number of optimum cell clusters for the given experimental scRNA-seq data.
#'
#' @param CountData Observed count data matrix for genes, rows represent genes, columns represent cells.
#' @param n Maximum value for number of cell clusters.
#' @param seed value for random cluster generation.
#' @param Threshold Threshold value for deciding the optimum number of cell clusters.
#' @param plot Logical variabe taking value either TRUE or FALSE, default is FALSE.
#'
#' @return A list containing the clustering index, delta and the optimum cluster number.
#'
#' @author Samarendra Das
#'
#' @examples
#' ##Load the test count data given for SwarnSeq.
#' counts <- matrix(rnbinom(2000, size=0.2, mu=3.4), 50)
#'
#' results <- optimcluster(CountData=counts, n = 10, seed = 108, Threshold = 0.3, plot = FALSE)
#'
#' @importFrom stats kmeans
#'
#' @export
#'
optimcluster <- function (CountData, n, seed, Threshold, plot=TRUE){
  set.seed(seed)
  if(is.null(n)) n <- 50
  if(is.null(plot)) plot = FALSE
  k <- seq(from=2, to=n, by=1)
  cl.ind <- vector(mode="numeric", length=length(k))
  for (i in 1:length(k)){
    cell_cl <- kmeans(t(CountData), centers=k[i])
    WSS <- cell_cl$tot.withinss
    BSS <- cell_cl$betweenss
    r <- WSS/BSS
    cl.ind[i] <- r
  }
  names(cl.ind) <- k
  delta <- vector()
  for (i in 1:length(cl.ind)) delta[i] = cl.ind[i] - cl.ind[i+1]
  delta <- delta[!is.na(delta)]; names(delta) <- k[-length(k)]
  Optim.cluster <- max(which(delta >  Threshold))
  if (plot==TRUE)
  {
    plot(k, cl.ind, typ="l", col="blue", lwd=2, ylab="Clustering Index", xlab="Number of clusters",
         xlim=c(2, max(k)), ylim=c(min(cl.ind), max(cl.ind)))
    #abline(v= Optim.cluster, lwd=2, col="red")
    #text(2, 1, "k=(Optim.cluster)", col = "red")
  }
  out <- list(cluster.index=cl.ind, delta=delta, Optimum.cluster=Optim.cluster)
  class(out) <- c("Clustering Index", "Optimum Cluster")
  return(out)
}

