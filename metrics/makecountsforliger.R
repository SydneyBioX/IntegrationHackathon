#Make sure my MOFA analysis R object is loaded into workspace
library(plyr)

genecountsfromatac <- function(gene){
  generange <- gtfranges[gene]
  ranges.a <- suppressWarnings(atac_ranges[atac_ranges %over% generange])
  if(length(ranges.a)==0){return(NULL)}
  if(length(ranges.a)==1){atacsums <- ataccounts[paste0(seqnames(ranges.a), "-", start(ranges.a), "-", end(ranges.a)),]
  } else {atacsums <- colSums(as.matrix(ataccounts[paste0(seqnames(ranges.a), "-", start(ranges.a), "-", end(ranges.a)),]))
  }
  atacsums
}

atac_liger <- t(sapply(unique(res$gene), genecountsfromatac))
expr_liger <- rnacounts[rownames(atac_liger), colnames(atac_liger)]