#' This function estimates the capture efficiencies of the cells from the single-cell RNA-seq studies. It takes input the ERCC spike-in transcript
#'  and molecular concentration data, if available or count expression data, if spike-ins are not available.

#' @param CountData Observed count data matrix for genes, rows represent genes, columns represent cells.
#' @param RNAspike.use Logical value indicating TRUE/FALSE, if TRUE, spikes and spike.conc information must be provided.
#' @param spikes Observed count matrix for spike-in transcripts, rows represent spike-in transcripts, columns represent cells. Only needed if RNAspike.use = \code{TRUE}).
#' @param spike.conc Vector of theoretical count for each spike-in transcript in one cell (ONLY needed if RNAspike.use = \code{TRUE}).
#' @param CE.range Two-element vector representing the lower and upper limits for the estimated range of capture efficiencies
#'  (ONLY needed if RNAspike.use = \code{FALSE}, default [0.1, 0.40]).
#' @param method Character representing the methods to be used for computation of capture efficiencies for cells.
#'
#' @return Returns a vector of estimated capture efficiencies for cells given in the scRNA-seq data.
#'
#' @author Samarendra Das
#'
#' @importFrom stats lm
#'
#' @export
#'
#'
CapEff <- function(CountData, CE.range, RNAspike.use, spikes, spike.conc, method){
  if (RNAspike.use) {
    if(method == "ML"){
      #spike1 <- sweep(spikes, 1, spike.conc, '*')
      capeff.spike <- apply(spikes, 2, sum)/sum(spike.conc)
      #capeff.spike <- apply(spike1, 2, sum)/sum(spike.conc)
      #DO.coef <- matrix(0, ncell, 2)
      #DO.coef[,1] <- log(capeff.spike/(1-capeff.spike))
      CE <- capeff.spike
      remove(spikes, capeff.spike);gc()
    }
    else{
      CE <- vector(mode="numeric", length=ncol(spikes))
      #R.CE <- vector(mode="numeric", length=ncol(spikes))
      for(i in 1:ncol(spikes)){
        spik <- as.numeric(spikes[,i])
        spike.conc <- as.vector(spike.conc)
        mod.spike <- lm(spik ~ spike.conc)
        CE[i] <- mod.spike$coefficients [2]
        names(CE) <- colnames(spikes)
        #R.CE[i] <- summary(mod.spike)$adj.r.squared
      }
      CE <- ifelse(CE < 0, 0, CE)
      if (CE >= 1) message("CE canot be more than 1; check inputs")
      remove(spik, spike.conc, mod.spike);gc()
    }
  } else {
    if(is.null(CE.range))
      CE.range <- c(0.01, 0.2)    #####user input
    if (CE.range[1] < 0 | CE.range[1] > CE.range[2] | CE.range[2] > 1) stop('CE.range invalid.')
    l.sz  <- log10(colSums(CountData))
    l.max  <- max(l.sz)
    l.min  <- min(l.sz)
    # generate rand.CE within CE.range but following the dist of obs.ls closely

    ls.wt <- (l.sz - l.min)/(l.max - l.min)
    rand.CE <- CE.range[1] + (CE.range[2] - CE.range[1]) * ls.wt
    #DO.coef <- matrix(0, ncell, 2)
    #DO.coef[, 1] <- log(CE/(1-CE))
    CE <- rand.CE; names(CE) <- colnames(CountData)

  }
  class(CE) <- "Capture Efficiency"
  return(CE)

}
