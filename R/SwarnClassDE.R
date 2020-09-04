#' This function is used to classify the influential genes of single-cell RNA-seq (scRNA-seq) data obtained from SwarnSeqLRT or SwarnUnadjLRT.
#'
#' @param results A output data frame from \code{SwarnSeqLRT} or \code{SwarnUnadjLRT} which contains the unclassified differential expression analysis results.
#' @param alpha A number in (0, 0.05) to specify the threshold of adjusted p-values.
#'
#' @return A list containing the results along with the classes of influential genes in scRNA-seq data.
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
#'
#' # Do not run
#' #results <- SwarnSeqLRT(CountData=counts, RNAspike.use=TRUE, spikes=Spikes, spike.conc=SpikeConc,
#'                       #parallel=FALSE, norm.method="TMM", group=group, CellCluster=cellcluster,
#'                           #CellAuxil=NULL, maxit=500, eps=1e-10,
#'                           #muoffset=NULL, phioffset=NULL, weights=NULL, p.adjust.method="hochberg")
#'
#' #DEGtypes <- SwarnClassDE(results, alpha = 0.0005)
#'
#' @export
#'
#'

SwarnClassDE <- function (results, alpha){

  # Invalid input judge
  if(class(results) != "data.frame")
    stop("Invalid input of wrong data type of results")
  if(ncol(results) != 16)
    stop("Invalid input of wrong column number of results, must be the object from 'SwarnSeqLRT'")
  if(length(results$DE.Adj.pval) != length(results$DZI.Adj.pval))
    stop("Invalid input of wrong column name of results")
  if(class(alpha) != "numeric")
    stop("Invalid input of wrong data type of alpha (Level of significance)")
  if(alpha <= 0 | alpha > 0.05)
    stop("Invalid input of wrong range of alpha (must be atleast 5%)")
  if(is.null(alpha))
    alpha <- 0.005
  # Classify the types of DE genes
  #results <- cbind(results, NA, NA)
  p.DE <- as.vector(results$DE.Adj.pval); names(p.DE) <- rownames(results) ###adjusted p values
  p.DZI <- as.vector(results$DZI.Adj.pval); names(p.DZI) <- rownames(results) # adjusted p values
  Class <- ifelse(p.DE < alpha & p.DZI < alpha, "DE&DZI",  ifelse(p.DE < alpha, "DE", ifelse(p.DZI < alpha, "DZI", "NonDE")))
  results <- cbind(results, Class)
  return(results)
}

###########
