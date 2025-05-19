## Lachlan sent me the code to get gene-level accessibility peak counts
## You would need to down load these two files from 10X (https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-3-k-1-standard-2-0-0)
## and save them under the same directory
## "pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
## "pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz.tbi"

library("ggplot2")
library("Seurat")
library("Signac")
library("MatrixGenerics")

fragments <- CreateFragmentObject(
  path = "data/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz", ## wherever it is
  cells = colnames(pbmc_multiome),
  validate.fragments = FALSE
)

# override the old one which has a path that does not exist
Fragments(pbmc_multiome) <- NULL
Fragments(pbmc_multiome) <- fragments

pbmc_multiome[["ATAC.GENE.ACTIVITY"]] <- CreateAssayObject(
  counts = GeneActivity(
    pbmc_multiome,
    features = rownames(pbmc_multiome@assays$RNA@layers$counts[rowSums(pbmc_multiome@assays$RNA@layers$counts) > 0,]))
)
