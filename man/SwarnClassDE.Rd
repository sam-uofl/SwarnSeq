% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SwarnClassDE.R
\name{SwarnClassDE}
\alias{SwarnClassDE}
\title{This function is used to classify the influential genes of single-cell RNA-seq (scRNA-seq) data obtained from SwarnSeqLRT or SwarnUnadjLRT.}
\usage{
SwarnClassDE(results, alpha)
}
\arguments{
\item{results}{A output data frame from \code{SwarnSeqLRT} or \code{SwarnUnadjLRT} which contains the unclassified differential expression analysis results.}

\item{alpha}{A number in (0, 0.05) to specify the threshold of adjusted p-values.}
}
\value{
A list containing the results along with the classes of influential genes in scRNA-seq data.
}
\description{
This function is used to classify the influential genes of single-cell RNA-seq (scRNA-seq) data obtained from SwarnSeqLRT or SwarnUnadjLRT.
}
\examples{

 #Load the test count data, spike-in counts and spike-in concentration data for SwarnSeq.
data(TestData)
counts <- TestData$CountData
Spikes <- TestData$SpikeCounts
SpikeConc <- TestData$SpikeConc

#specifying the group information, the group 1 and 2 have two hundred cells each.
group <- c(rep(1, 200), rep(2, 200))
#Specifying the cluster memberships of the cells in columns of countData.
cellcluster <- c(rep(1, 60), rep(2, 40), rep(3, 50),
                    rep(4, 50), rep(5, 30),
                    rep(6, 90),
                    rep(7, 80))

#results <- SwarnSeqLRT(CountData=counts, RNAspike.use=TRUE, spikes=Spikes, spike.conc=SpikeConc,
                      #parallel=FALSE, norm.method="TMM", group=group, CellCluster=cellcluster,
                          #CellAuxil=NULL, maxit=500, eps=1e-10,
                          #muoffset=NULL, phioffset=NULL, weights=NULL, p.adjust.method="hochberg")

#DEGtypes <- SwarnClassDE(results, alpha = 0.0005)

}
\seealso{
\code{\link{SwarnSeqLRT}}, for the detection of differentially expressed genes from scRNA-seq data.

\code{\link{TestData}}, a test dataset for SwarnSeq.
}
\author{
Samarendra Das
}
