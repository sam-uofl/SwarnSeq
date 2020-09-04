#' This function is used to detect differentially expressed genes between two specified groups of cells in a raw read counts matrix of single-cell RNA-seq (scRNA-seq) data without adjustment for capture efficiency. It takes a non-negative integer matrix of scRNA-seq raw read counts object as input. So users should map the reads (obtained from sequencing libraries of the samples) to the corresponding genome and count the reads mapped to each gene according to the gene annotation to get the raw read counts matrix in advance.
#'
#' @param CountData Observed count data matrix for genes, rows represent genes, columns represent cells.
#' @param parallel If FALSE (default), no parallel computation is used; if TRUE, parallel computation is performed.
#' @param norm.method Method for normalizing the scRNA-seq count expression data, either 'DESeq2' (maximum likelihood, Ye et al., 2017) or 'TMM' (Robinson et al., 2010).
#' @param group Vector which specifies the membership of the cells, i.e. two groups to be compared, corresponding to the columns in the count data matrix.
#' @param CellCluster Vector which specifies the cluster memberships of the cells, i.e. each entry represents memberships of the columns of the count data matrix.
#' @param CellAuxil Vector of cell level auxiliary information, corresponding to the columns in the counts matrix, default is NULL.
#' @param maxit Maximum number of iterations for Expected-Maximization (EM) algorithm.
#' @param eps Convergennce criteria for EM algorithm.
#' @param muoffset Offset parameter for mean (mu) parameter, default is NULL.
#' @param phioffset Offset parameter for zero inflation (phi) parameter, default is NULL.
#' @param weights Observation wise weights for the cells, default is unity vector.
#' @param p.adjust.method Character variable represents the method used for multiple hypothesis correction. It can be any value from ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY").
#'
#' @return
#' A data frame containing the results from differential expression analysis, rows are genes and columns contain the following items:
#' \itemize{
#'  \item1 totalMean_1, totalMean_2,	NormMean_1,	NormMean_2 are the total mean, normalized mean for cellular groups 1 and 2 respectively.
#'  \item2 FoldChange,	log2FC,	NormFC,	log2NormFC are the fold change, log fold change, and log normalized fold change for the genes respectively.
#'  \item3 Stat.DE,	Pval.DE,	DE.Adj.pval, and	DE.FDR are values of DE statistic, p-value, adjusted p-value, false discovery rate, obtained from DE analysis, for the genes.
#'  \item4 Stat.DZI,	Pval.DZI,	DZI.Adj.pval,	DZI.FDR are Differential Zero Inflation (DZI) statistic, DZI p-value, DZI adjusted p-value, DZI false discovery rate results obtained for each gene from DZI analysis.
#'  }
#'
#' @author Samarendra Das
#'
#' @seealso \code{\link{TestData}}, a test dataset for SwarnSeq.
#' @seealso \code{\link{dzinb}}, dzinb funcion in SwarnSeq.
#' @seealso \code{\link{ZINBEM}}, Expected-Maximization (EM) algorithm in SwarnSeq.
#' @seealso \code{\link{SwarnUnadjSeq}}, SwarnUnadjSeq function in SwarnSeq.
#'
#' @examples
#'
#' #Load the test count data, spike-in counts and spike-in concentration data for SwarnSeq.
#' data(TestData)
#' counts <- TestData$CountData
#'
#' #specifying the group information, the group 1 and 2 have two hundred cells each.
#' group <- c(rep(1, 200), rep(2, 200))
#' #Specifying the cluster memberships of the cells in columns of countData.
#' cellcluster <- c(rep(1, 60), rep(2, 40), rep(3, 50),
#'                     rep(4, 50), rep(5, 30),
#'                     rep(6, 90),
#'                     rep(7, 80))
#'
#' #Do not run
#' #parameters from EM algorithm for each gene.
#' #results <- SwarnUnadjLRT(CountData=counts, parallel=FALSE, norm.method="TMM", group=group,
#'                       # CellCluster=cellcluster, CellAuxil=NULL, maxit=500, eps=1e-10,
#'                           #muoffset=NULL, phioffset=NULL, weights=NULL, p.adjust.method="hochberg")
#'
#' @import stats
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb
#' @importMethodsFrom Matrix colSums
#' @importFrom foreach foreach %do%
#' @importFrom edgeR calcNormFactors
#'
#' @export
#'
#'
SwarnUnadjLRT <- function(CountData, parallel, norm.method, group, CellCluster, CellAuxil, maxit, eps, muoffset, phioffset, weights, p.adjust.method){

  if(!is.matrix(CountData) & !is.data.frame(CountData) & class(CountData)[1] != "dgCMatrix")
    stop("Wrong input data type of 'Count Data'")
  if(sum(is.na(CountData)) > 0)
    stop("NAs detected in input 'Count Data'");gc();
  if(sum(CountData < 0) > 0)
    stop("Negative values detected in input 'Count Data'");gc();
  if(all(CountData == 0))
    stop("All elements of input 'Count Data' are zeros");gc();
  if(any(colSums(CountData) == 0))
    warning("Library size of zeros detected in 'Counts Data'");gc();
  #if(!is.factor(group) )
  #stop("Data type of 'group' is not factor")
  if(max(group) != 2)
    stop("Levels number of 'group' is not two")
  if(table(group)[1] < 2 | table(group)[2] < 2)
    stop("Too few samples (< 2) in a group")
  if(ncol(CountData) != length(group) | ncol(CountData) != length(CellCluster))
    stop("Length of 'group' and 'CellCluster' must equal to column number of 'Count Data'")
  if(!is.null(CellAuxil) & ncol(CountData) != length(CellAuxil))
    stop ("Length of cell_level co-variate factor should be same as number of cells")

  if(!is.logical(parallel))
    stop("Data type of 'parallel' is not logical")
  if(length(parallel) != 1)
    stop("Length of 'parallel' is not one")
  if(!is.numeric(c(maxit, eps)))
    stop("Data type of 'maxit' & 'eps' is not numeric")
  if(length(maxit) != 1)
    stop("Length of 'maxit' is not one")
  if(length(eps) != 1)
    stop("Length of converegence criterion 'eps' is not one")
  if(!is.character(norm.method))
    stop("Data type of 'Norm.method' is not character")
  if(!is.element(norm.method, c("DEseq.norm", "TMM")))
    stop("Normalization method mlust be either from  'DEseq.norm', 'TMM'")

  if(!is.character(p.adjust.method))
    stop("Data type of 'p.adjust.method' is not character")
  method <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY")
  if(!is.element(p.adjust.method, method))
    stop("Multiple hypothesis correction method must be either from  'holm', 'hochberg',
         'hommel', 'bonferroni', 'BH', 'BY'")
  if(is.null(p.adjust.method))
    p.adjust.method <- "BH"
  #######

  ####pre-processing
  CountData <- round(as.matrix(CountData))
  storage.mode(CountData) <- "integer"
  if(any(rowSums(CountData) == 0))
    message("Removing ", sum(rowSums(CountData) == 0), "Genes with zero counts across all the cells removed...")
  CountData <- CountData[rowSums(CountData)!=0,]

  f0 <- function(y) length(y) - sum(y==0)
  fg <- apply(CountData, 1, f0)
  if(any(fg < 5))
    message("Removing ", sum(fg <= 5),"\t","Genes, which dont have atleat five non-zero counts across all the cells ...")
  CountData <- CountData[fg > 5, ]
  Ngene <- nrow(CountData)
  Ncell <- ncol(CountData)
  remove(fg); gc()

  group <- as.factor(group)
  ############total mean and fold chage

  # Cache totalMean and foldChange for each gene
  totalMean_1 <- rowMeans(CountData[row.names(CountData), group == levels(group)[1]])
  totalMean_1[totalMean_1==0] <- 0.1
  totalMean_2 <- rowMeans(CountData[row.names(CountData), group == levels(group)[2]])
  totalMean_2[totalMean_2==0] <- 0.1
  foldChange <- totalMean_2 / totalMean_1
  logFC <- log2(foldChange)
  All_Mean_FC <- cbind(totalMean_1, totalMean_2, foldChange, logFC); colnames(All_Mean_FC) <- c("TotalMean_1","TotalMean_2","FoldChange", "log2FC")
  ######memeory management
  remove(totalMean_1, totalMean_2, foldChange, logFC)
  #########Normalization of the count data
  cat("\n", paste("Normalizing the scRNA-seq read counts..."))
  #NormalizeCountData <- function(countData, norm.method)
  #{
  if (norm.method=="DEseq.norm")
  {
    GM <- function(x) exp(mean(log(x[x>0]))) ######geometric mean
    geomMean <- apply(CountData, 1, GM)
    fx <- function(x) x/geomMean
    samp <- apply (CountData, 2, fx)
    fxx <- function(xx) median(xx[xx!=0])
    size <- apply(samp, 2, fxx)
    fxxx <- function(xxx) xxx/size
    counts_norm <- t(apply(CountData, 1, fxxx))
    counts_norm <- ceiling(counts_norm)
    #return(counts_norm)
    remove(CountData, geomMean, samp, size)
    gc()
  }
  if(norm.method=="TMM")
  {
    size_factor <- suppressWarnings(edgeR::calcNormFactors(CountData))
    fxxx <- function(xxx) xxx/size_factor
    counts_norm <- t(apply(CountData, 1, fxxx))
    counts_norm <- ceiling(counts_norm)
    #return(counts_norm)
    remove(size_factor, CountData);gc()
  }
  #else if(normalization=="DESeq2_poscounts"){
  #dse = DESeqDataSetFromMatrix(countData, colData=colData, design=designFormula)
  #dse = estimateSizeFactors(dse, type = "poscounts")
  #counts$samples$norm.factors = 1/dse$sizeFactor

  # Cache totalMean and foldChange for each gene
  totNormMean_1 <- rowMeans(counts_norm[row.names(counts_norm), group == levels(group)[1]])
  totNormMean_1[totNormMean_1 == 0] <- 0.1
  totNormMean_2 <- rowMeans(counts_norm[row.names(counts_norm), group == levels(group)[2]])
  totNormMean_2[totNormMean_2 == 0] <- 0.1
  NormFC <- totNormMean_2/totNormMean_1
  logNormFC <- log2(NormFC)
  Norm_Mean_FC <- cbind(totNormMean_1, totNormMean_2, NormFC, logNormFC); colnames(Norm_Mean_FC) <- c("NormMean_1","NormMean_2","NormFC", "log2NormFC")
  # Memory management
  remove(totNormMean_1, totNormMean_2, NormFC, logNormFC)
  counts_norm <- Matrix(counts_norm, sparse = TRUE)
  gc()
  #}
  #counts_norm <- CountData
  # Call DEG gene by gene
  if(!parallel){
    results <- matrix(data=NA, nrow = Ngene, ncol = 4)
    results <- as.data.frame(results)
    for(i in 1: Ngene){
      cat("\r",paste0("SwarnSeq is analyzing ", i," of ",Ngene," expressed genes..."))
      Y <- as.numeric(counts_norm[i,])
      EM_Gene <- ZINBEM (Count=Y, group, CellCluster, CellAuxil, maxit, eps, muoffset, phioffset, weights)
      results[i,] <- EM_Gene$DE.stat; colnames(results) <- names(EM_Gene$DE.stat)
    }
  }else{
    callDE <- function(i){
      reslt <- ZINBEM (Count=as.numeric(counts_norm[i,]), group, CellCluster, CellAuxil, maxit, eps, muoffset, phioffset, weights)$DE.stat
      reslt
    }
    cat("\n", paste("SwarnSeq is analyzing ", Ngene, " expressed genes in parallel..."))
    results <- foreach(i=1:Ngene, .combine=rbind) %do% callDE(i)
    rownames(results) <- rownames(counts_norm)
  }

  ########Adjusting p-values
  cat("\n", paste("Correction for multiple hypothesis testing..."))
  adj.method <- method[match(p.adjust.method, method)]
  if(is.null(p.adjust.method))
    adj.method <- "BH"
  p.DE <- as.vector(results[, 2])
  adj.pval.DE <- p.adjust(p.DE, method = adj.method, n = length(p.DE))
  FDR.DE <-  p.adjust(p.DE, method = "fdr", n = length(p.DE))
  p.DZI <- as.vector(results[, 4])
  adj.pval.DZI <- p.adjust(p.DZI, method = adj.method, n = length(p.DZI))
  FDR.DZI <-  p.adjust(p.DZI, method = "fdr", n = length(p.DZI))
  p.adj <- cbind(adj.pval.DE, FDR.DE, adj.pval.DZI, FDR.DZI); colnames(p.adj) <- c("DE.Adj.pval", "DE.FDR", "DZI.Adj.pval", "DZI.FDR")
  remove(p.DE, p.DZI, adj.pval.DE, FDR.DE, adj.pval.DZI, FDR.DZI)
  results <- cbind(All_Mean_FC,  Norm_Mean_FC, results, p.adj)
  remove(All_Mean_FC, Norm_Mean_FC, p.adj); gc()
  return(results)
}
