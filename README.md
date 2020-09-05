# SwarnSeq
This R package performs Differential Expression, and Differential Zero Inflation analysis of  Single-Cell RNA-seq (scRNA-seq) data.

Differential Expression (DE) analysis is one of the powerful downstream analysis for scRNA-seq data.
It is also required for obtaining informative gene sets and further used as input for other analysis such as pathway or gene set analysis, 
gene regulation modeling, gene network analysis, etc. Further, there is a lot of methods and tools are available in literature for performing DE analysis of scRNA-seq data, 
which significantly differ from each other on the basis of statistical methods, input data, etc. The existing methods mostly ignore the cell level auxiliary information 
such as cell cluster, cell cycle, molecular capturing procedure, etc. and limited to only two cellular populations comparison from the model building. 
For instance, DEsingle, a popular method does not consider molecular capturing process, cellular cluster, cell cycle, cell phase, etc. for DE analysis and only limited
to two cell groups comparisons. Therefore, proposed SwarnSeq approach and the R package performs various analysis, such as DE analysis, Differential Zero Inflation analysis and classification of genes
into subtypes, of scRNA-seq data by considering RNA capturing process, cell level auxiliary information (i.e. cell cluster, other cell level data) in the statistical model building. 
