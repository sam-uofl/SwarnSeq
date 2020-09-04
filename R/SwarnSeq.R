#'
#' This function is used to estimate the parameters for genes between two specified groups of cells in a raw read counts matrix of single-cell RNA-seq (scRNA-seq) data. It takes a non-negative integer matrix of scRNA-seq raw read counts as input. So users should map the reads (obtained from sequencing libraries of the samples) to the corresponding genome and count the reads mapped to each gene according to the gene annotation to get the raw read counts matrix in advance.
#'
#'
#' @param CountData Observed count data matrix for genes, rows represent genes, columns represent cells.
#' @param RNAspike.use Logical value indicating TRUE/FALSE, if TRUE, spikes and spike.conc information must be provided.
#' @param spikes Observed count matrix for spike-in transcripts, rows represent spike-in transcripts, columns represent cells. Only needed if RNAspike.use = \code{TRUE}).
#' @param spike.conc Vector of theoretical count for each spike-in transcript in one cell (ONLY needed if RNAspike.use = \code{TRUE}).
#' @param CE.range Two-element vector representing the lower and upper limits for the estimated range of capture efficiencies
#'  (ONLY needed if RNAspike.use = \code{FALSE}, default [0.1, 0.40]).
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
#'
#' @return
#' A data frame containing the parameters from the EM algorithm for each gene, rows are genes and columns contain the following items:
#'\itemize{
#'  \item1 totalMean_1, totalMean_2,	Adj_Mean_1,	Adj_Mean_2,	AdjNormMean_1,	AdjNormMean_2 are the total mean, adjusted mean, normalized mean for cellular groups 1 and 2 respectively.
#'  \item2 Mean, ZeroInflation,	Dispersion are the MLE of the parameters of whole cellular population characetrized by the ZINB model.
#'  \item3 Intercept,	Group2 are the co-efficients of the intercept and group effect on mean of non-zero counts.
#'  \item4 Cellcluster 2, Cellcluster 3, .., Cellcluster M are the effects of cell clusters on mean of non-zero counts.
#'  \item5 CellAuxil 2, CellAuxil 3, ..., CellAuxil N are the effects of cell-level auxiliary information on mean of non-zero counts (if included in the model).
#'  \item6 Intercept.0	Group2.0 are the co-efficients of the intercept and group effect on zero inflation probability.
#'  \item7 Cellcluster 2.0, Cellcluster 3.0, .., Cellcluster M.0 are the effects of cell clusters on zero-inflation probability.
#'  \item8 CellAuxil 2.0, CellAuxil 3.0, ..., CellAuxil N.0 are the effects of cell-level auxiliary information on zero-inflation probability (if included in the model).
#'  \item9 #iteration number of iteration required for convergence for each gene.
#'  }
#'
#' @author Samarendra Das
#'
#' @seealso \code{\link{TestData}}, a test dataset for SwarnSeq.
#' @seealso \code{\link{dzinb}}, dzinb funcion in SwarnSeq.
#' @seealso \code{\link{ZINBEM}}, Expected-Maximization (EM) algorithm in SwarnSeq.
#'
#'
#' @examples
#'
#' #Load the test count data, spike-in counts and spike-in concentration data for SwarnSeq.
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
#' #parameters from EM algorithm for each gene (when Spike-in information available).
#' #results <- SwarnSeq(CountData=counts, RNAspike.use=TRUE, spikes=Spikes, spike.conc=SpikeConc,
#'                       #parallel=FALSE, norm.method="TMM", group=group, CellCluster=cellcluster,
#'                           #CellAuxil=NULL, maxit=500, eps=1e-10,
#'                           #muoffset=NULL, phioffset=NULL, weights=NULL)
#'
#' #When Spike-in information not available.
#' #results <- SwarnSeq(CountData=counts, RNAspike.use=FALSE, CE.range=c(0.1, 0.4)
#'                       #parallel=FALSE, norm.method="TMM", group=group, CellCluster=cellcluster,
#'                           #CellAuxil=NULL, maxit=500, eps=1e-10,
#'                           #muoffset=NULL, phioffset=NULL, weights=NULL)
#'
#' @importFrom stats lm model.matrix median glm.fit optim
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb
#' @importMethodsFrom Matrix colSums
#' @importFrom foreach foreach %do%
#' @importFrom edgeR calcNormFactors
#'
#' @export
#'
#'
#'
SwarnSeq <- function(CountData, RNAspike.use, spikes, spike.conc, CE.range, parallel, norm.method, group, CellCluster, CellAuxil, maxit, eps, muoffset, phioffset, weights){
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

  ####Optimization, when EM algorthm failed

  ##############EM algorithm

  ###########co-variates
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
  kx <- NCOL(X)
  kz <- NCOL(Z)
  remove(X, Z); gc()

  ####pre-processing
  CountData <- round(as.matrix(CountData))
  storage.mode(CountData) <- "integer"
  if(any(rowSums(CountData) == 0))
    message("Removing ", sum(rowSums(CountData) == 0),"\t","Genes with zero counts across all the cells removed...")
  CountData <- CountData[rowSums(CountData)!=0,]

  f0 <- function(y) length(y) - sum(y==0)
  fg <- apply(CountData, 1, f0)
  if(any(fg < 5))
    message("Removing ", sum(fg <= 5),"\t","Genes, which dont have atleat five non-zero counts across all the cells ...")
  CountData <- CountData[fg > 5, ]
  Ngene <- nrow(CountData)
  Ncell <- ncol(CountData)
  remove(fg); gc()
  ############Estimation of capture efficiency
  #cell.type <- rep(1, ncell)
  #cell.type.names <- NULL
  #ncelltype <- length(unique(cell.type))
  #XW <- as.matrix(model.matrix(~cell.type)[,1])
  # Get capture efficiency. Calculate with spike-ins, if available;
  # If not, randomize and sort by libsize.
  ######Estimation of capture efficiency
  cat (paste ("Adjusting scRNA-seq reads for capture efficiency..."))
  if (RNAspike.use) {
    CE <- vector(mode="numeric", length=ncol(spikes))
    for(i in 1:ncol(spikes)){
      spik <- as.numeric(spikes[,i])
      spike.conc <- as.numeric(spike.conc)
      mod.spike <- lm(spik ~ spike.conc)
      CE[i] <- mod.spike$coefficients [2]
    }
    CE <- ifelse(CE<0, 0, CE)
    remove(spikes, spik, spike.conc, mod.spike);gc()
  } else {
    l.sz  <- log10(colSums(CountData))
    l.max  <- max(l.sz)
    l.min  <- min(l.sz)
    # generate rand.CE within CE.range but following the dist of obs.ls closely
    if(is.null(CE.range))
      CE.range <- c(0.01, 0.2)    #####user input
    ls.wt <- (l.sz - l.min)/(l.max - l.min)
    rand.CE <- CE.range[1] + (CE.range[2] - CE.range[1]) * ls.wt
    #DO.coef <- matrix(0, ncell, 2)
    #DO.coef[, 1] <- log(CE/(1-CE))
    CE <- as.vector(rand.CE)
    remove(ls.wt, rand.CE, l.sz)
  }
  counts.adj <- sweep(CountData, 2, CE, '/')
  ####Alternatively
  #counts.adj = CountData %*% diag(1/CE)
  #remove(ls.wt, rand.CE, l.sz)

  # Cache totalMean and foldChange for each gene
  totalMean_1 <- rowMeans(CountData[row.names(CountData), group == levels(group)[1]])
  totalMean_1[totalMean_1==0] <- 0.1
  totalMean_2 <- rowMeans(CountData[row.names(CountData), group == levels(group)[2]])
  totalMean_2[totalMean_2==0] <- 0.1
  All_Mean <- cbind(totalMean_1, totalMean_2); colnames(All_Mean) <- c("totalMean_1", "totalMean_2")
  ######memeory management
  remove(CountData, totalMean_1, totalMean_2); gc()

  #############cache Total mean for geans for adjusted scRNA-seq data
  Adj_Mean_1 <- rowMeans(counts.adj[row.names(counts.adj), group == levels(group)[1]])
  Adj_Mean_1[Adj_Mean_1 == 0] <- 0.1
  Adj_Mean_2 <- rowMeans(counts.adj[row.names(counts.adj), group == levels(group)[2]])
  Adj_Mean_2[Adj_Mean_2 == 0] <- 0.1
  All_Adj_Mean <- cbind(Adj_Mean_1, Adj_Mean_2); colnames(All_Adj_Mean) <- c("Adj_Mean_1", "Adj_Mean_2")
  #####Memory management
  remove(Adj_Mean_1, Adj_Mean_2); gc()
  #########Normalization of the count data


  #NormalizeCountData <- function(countData, norm.method)
  #{
  if (norm.method=="DEseq.norm")
    cat(paste("\n", "Normalizing the adjusted scRNA-seq read counts using DEseq.norm..."))
  {
    GM <- function(x) exp(mean(log(x[x>0]))) ######geometric mean
    geomMean <- apply(counts.adj, 1, GM)
    fx <- function(x) x/geomMean
    samp <- apply (counts.adj, 2, fx)
    fxx <- function(xx) median(xx[xx!=0])
    size <- apply(samp, 2, fxx)
    fxxx <- function(xxx) xxx/size
    counts_norm <- t(apply(counts.adj, 1, fxxx))
    counts_norm <- ceiling(counts_norm)
    #return(counts_norm)
    remove(geomMean, samp, size)
    gc()
    }
  if(norm.method=="TMM")
  {
    cat(paste("\n", "Normalizing the adjusted scRNA-seq read counts using TMM..."))
    size_factor <- suppressWarnings(calcNormFactors(counts.adj))
    fxxx <- function(xxx) xxx/size_factor
    counts_norm <- t(apply(counts.adj, 1, fxxx))
    counts_norm <- ceiling(counts_norm)
    #return(counts_norm)
    remove(size_factor, counts.adj);gc()
  }
  #else if(normalization=="DESeq2_poscounts"){
  #dse = DESeqDataSetFromMatrix(countData, colData=colData, design=designFormula)
  #dse = estimateSizeFactors(dse, type = "poscounts")
  #counts$samples$norm.factors = 1/dse$sizeFactor

  # Cache totalMean and foldChange for each gene
  totAdjNormMean_1 <- rowMeans(counts_norm[row.names(counts_norm), group == levels(group)[1]])
  totAdjNormMean_1[totAdjNormMean_1 == 0] <- 0.1
  totAdjNormMean_2 <- rowMeans(counts_norm[row.names(counts_norm), group == levels(group)[2]])
  totAdjNormMean_2[totAdjNormMean_2 == 0] <- 0.1
  Adj_Norm_Mean <- cbind(totAdjNormMean_1, totAdjNormMean_2); colnames(Adj_Norm_Mean) <- c("AdjNormMean_1", "AdjNormMean_2")
  # Memory management
  remove(totAdjNormMean_1, totAdjNormMean_2); gc()
  counts_norm <- Matrix(counts_norm, sparse = TRUE)
  gc()
  #}

  # Call DEG gene by gene
  if(!parallel){
    results <- matrix(data=NA, nrow = Ngene, ncol = (3 + kx + kz + 1))
    results <- as.data.frame(results)
    for(i in 1: Ngene){
      cat("\n",paste("SwarnSeq is analyzing ", i," of ",Ngene," expressed genes ..."))
      Y <- as.numeric(counts_norm[i,])
      EM_Gene <- ZINBEM (Count=Y, group = group, CellCluster = CellCluster, CellAuxil = CellAuxil, maxit = maxit, eps = eps, muoffset = muoffset,
                          phioffset = phioffset, weights = weights)
      results[i,] <- EM_Gene$res; colnames(results) <- names(EM_Gene$res)
      rownames(results) <- rownames(counts_norm)
      #######memory management
      remove(Y); gc()
    }
  }else{
    callDE <- function(i){
      reslt <- ZINBEM(Count=as.numeric(counts_norm[i,]), group = group, CellCluster = CellCluster, CellAuxil = CellAuxil, maxit = maxit, eps = eps,
                       muoffset = muoffset, phioffset = phioffset, weights= weights)$res
      #reslt <- EM_Gene$res; names(reslt) <- names(EM_Gene$res)
      reslt
    }
    cat("\n", paste("SwarnSeq is analyzing ", Ngene, " expressed genes in parallel..."))
    results <- foreach(i=1:Ngene, .combine=rbind) %do% callDE(i)
    rownames(results) <- rownames(counts_norm)
  }
  results <- cbind(All_Mean, All_Adj_Mean, Adj_Norm_Mean, results)
  #######memory management
  remove(All_Mean, All_Adj_Mean, Adj_Norm_Mean); gc()
  return(results)
}

#################
