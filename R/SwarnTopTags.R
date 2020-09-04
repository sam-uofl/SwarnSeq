#' This function selects the top genes in scRNA-seq data.
#'
#' @param results A output data frame from \code{SwarnSeqLRT} or \code{SwarnUnadjLRT} which contains the unclassified differential expression analysis results.
#' @param m A scalar representing the number of top performing genes to be selected from the scRNA-seq data.
#'
#' @return A list of the top genes along with their statistics.
#'
#' @author Samarendra Das
#'
#' @seealso
#' \code{\link{SwarnSeqLRT}}, for the detection of differentially expressed genes from scRNA-seq data.
#'
#' \code{\link{TestData}}, a test dataset for SwarnSeq.
#'
#' @examples
#'
#'  #Load the test count data, spike-in counts and spike-in concentration data for SwarnSeq.
#' data(TestData)
#' counts <- TestData$CountData
#' Spikes <- TestData$SpikeCounts
#' SpikeConc <- TestData$SpikeConc
#'
#' #specifying the group information, the group 1 and 2 have two hundred cells each.
#' group <- c(rep(1, 200), rep(2, 200))
#' #Specifying the cluster memberships of the cells in columns of countData.
#' cellcluster <- c(rep(1, 60), rep(2, 40), rep(3, 50),
#'                     rep(4, 50), rep(5, 30),
#'                     rep(6, 90),
#'                     rep(7, 80))
#' #Do not run
#' #results <- SwarnSeqLRT(CountData=counts, RNAspike.use=TRUE, spikes=Spikes, spike.conc=SpikeConc,
#'                       #parallel=FALSE, norm.method="TMM", group=group, CellCluster=cellcluster,
#'                           #CellAuxil=NULL, maxit=500, eps=1e-10,
#'                           #muoffset=NULL, phioffset=NULL, weights=NULL, p.adjust.method="hochberg")
#'
#' #TopGene <- SwarnTopTags(results, m = 100)
#'
#' @export
#'

SwarnTopTags <- function (results, m){

  # Invalid input judge
  if(class(results) != "data.frame")
    stop("Invalid input of wrong data type of results")
  if(ncol(results) != 16)
    stop("Invalid input of wrong column number of results, must be the object from 'SwarnSeqLRT'")
  if(length(results$DE.Adj.pval) != length(results$DZI.Adj.pval))
    stop("Invalid input of wrong column name of results")
  if(class(m) != "numeric")
    stop("Invalid input of wrong data type of n (number of top tags)")
  if(m <= 0 | m > nrow(results))
    stop("Invalid input of wrong value of m")

  #####Top Tags

  p.DE <- as.vector(results$DE.Adj.pval); names(p.DE) <- rownames(results) ###adjusted p values
  p.DZI <- as.vector(results$DZI.Adj.pval); names(p.DZI) <- rownames(results) # adjusted p values
  id <- sort(p.DE, decreasing = FALSE, index.return = TRUE)$ix
  id <- id[1:m]
  top.DE <- results[id, ]
  id1 <- sort(p.DZI, decreasing = FALSE, index.return = TRUE)$ix
  id1 <- id1[1:m]
  top.DZI <- results[id1, ]
  out <- list(Top.DEG = top.DE, Top.DZIG = top.DZI)
  remove(results, id1, id, p.DE, p.DZI); gc()
  class(out) <- c("Top DE and zero inflated Tags")
  return(out)
}

#########
