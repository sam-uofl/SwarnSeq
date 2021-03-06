% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimcluster.R
\name{optimcluster}
\alias{optimcluster}
\title{This function decides the number of optimum cell clusters for the given experimental scRNA-seq data.}
\usage{
optimcluster(CountData, n, seed, Threshold, plot = TRUE)
}
\arguments{
\item{CountData}{Observed count data matrix for genes, rows represent genes, columns represent cells.}

\item{n}{Maximum value for number of cell clusters.}

\item{seed}{value for random cluster generation.}

\item{Threshold}{Threshold value for deciding the optimum number of cell clusters.}

\item{plot}{Logical variabe taking value either TRUE or FALSE, default is FALSE.}
}
\value{
A list containing the clustering index, delta and the optimum cluster number.
}
\description{
This function decides the number of optimum cell clusters for the given experimental scRNA-seq data.
}
\examples{
##Load the test count data given for SwarnSeq.
counts <- matrix(rnbinom(2000, size=0.2, mu=3.4), 50)

results <- optimcluster(CountData=counts, n = 10, seed = 108, Threshold = 0.3, plot = FALSE)

}
\author{
Samarendra Das
}
