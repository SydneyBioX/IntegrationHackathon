library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(scuttle)
library(scran)

sobj <- readRDS("~/analysis_tmp/Peters_workshop/pbmc_multiome.RDS")

DefaultAssay(sobj) <- "RNA"
#sobj <- NormalizeData(object = sobj)
#sobj <- FindVariableFeatures(object = sobj)
#sobj <- ScaleData(object = sobj)
sobj <- SCTransform(sobj)
sobj <- RunPCA(object = sobj)
#sobj <- FindNeighbors(object = sobj, dims = 1:30)
#sobj <- FindClusters(object = sobj)
sobj <- RunUMAP(object = sobj, dims = 1:20)

#Annotate with Monaco reference
rnacounts <- as.matrix(sobj@assays$RNA$counts)
rnacounts <- rnacounts[rowSums(rnacounts) > 0,]
sce <- SingleCellExperiment(list(counts = rnacounts), 
                            rowData = rownames(rnacounts), 
                            colData = sobj@meta.data)
sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)
logCPM <- logcounts(sce)

logCPM_RNA <- logCPM

#Annotate with Monaco 

library(SingleR)
library(celldex)
library(ggplot2)
ref <- MonacoImmuneData()

pred.cells <- SingleR(test = sce, ref = ref, labels = ref$label.main)

sobj <- AddMetaData(sobj, pred.cells$labels, col.name="Celltype")
DimPlot(sobj, reduction = "umap", group.by = "Celltype") + ggtitle("RNA UMAP")

#ATAC UMAP
DefaultAssay(sobj) <- "ATAC"
#sobj <- NormalizeData(object = sobj)
#sobj <- FindVariableFeatures(object = sobj)
#sobj <- ScaleData(object = sobj)
sobj <- SCTransform(sobj, assay = "ATAC")
sobj <- RunPCA(object = sobj)
#sobj <- FindNeighbors(object = sobj, dims = 1:30)
#sobj <- FindClusters(object = sobj)
sobj <- RunUMAP(object = sobj, dims = 1:20)
DimPlot(sobj, reduction = "umap", group.by = "Celltype") + ggtitle("ATAC UMAP")

#Back to RNA

DefaultAssay(sobj) <- "RNA"
#sobj <- NormalizeData(object = sobj)
#sobj <- FindVariableFeatures(object = sobj)
#sobj <- ScaleData(object = sobj)
sobj <- SCTransform(sobj)
sobj <- RunPCA(object = sobj)
#sobj <- FindNeighbors(object = sobj, dims = 1:30)
#sobj <- FindClusters(object = sobj)
sobj <- RunUMAP(object = sobj, dims = 1:20)

pdf("~/analysis_tmp/Peters_workshop/RNA_top100_mvgs.pdf", 12, 8)
for (i in 1:25){
  genes <- sobj@assays$SCT@var.features[((i-1)*4+1):(i*4)]
  p1 <- as.grob(FeaturePlot(sobj, features=genes[1], pt.size = 0.2) + ggtitle(genes[1]))
  p2 <- as.grob(FeaturePlot(sobj, features=genes[2], pt.size = 0.2) + ggtitle(genes[2]))
  p3 <- as.grob(FeaturePlot(sobj, features=genes[3], pt.size = 0.2) + ggtitle(genes[3]))
  p4 <- as.grob(FeaturePlot(sobj, features=genes[4], pt.size = 0.2) + ggtitle(genes[4]))
  grid.arrange(grobs=list(p1, p2, p3, p4), cols=2)
  print(i)
}
dev.off()

ataccounts <- sobj@assays$ATAC$counts
sce <- SingleCellExperiment(list(counts = ataccounts), 
                              rowData = rownames(ataccounts), 
                              colData = sobj@meta.data)
sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)
logCPM <- logcounts(sce)
  
logCPM_ATAC <- logCPM

##############MOFA##################

#At least 5% expressed or open
expr <- logCPM_RNA[apply(logCPM_RNA, 1, function (x) (sum(x > 0)/length(x)) > 0.05),]

atac <- logCPM_ATAC[apply(logCPM_ATAC, 1, function (x) (sum(x > 0)/length(x)) > 0.05),]
atac <- as.matrix(atac)

library(reticulate)
py_discover_config()
#virtualenv_install(~/.virtualenvs/r-reticulate/, python=/usr/bin/python3.7, packages=mofapy2)
#use_python(/usr/bin/python3.7)
library(MOFA2)
library(MultiAssayExperiment)


ExpList <- ExperimentList()
assayList <- list(EXPR=expr, ATAC=atac)
ExpList <- ExperimentList(assayList)
exprmap <- data.frame(primary=colnames(expr), colname=colnames(expr), stringsAsFactors = F)
atacmap <- data.frame(primary=colnames(atac), colname=colnames(atac), stringsAsFactors = F)
maplist <- list(EXPR=exprmap, ATAC=atacmap)
sampmap <- listToMap(maplist)
pd <- sobj@meta.data
pd$primary <- colnames(expr)
rownames(pd) <- pd$primary

GC_MOFA <- MultiAssayExperiment(experiments = ExpList, colData = pd, sampleMap = sampmap)
GC_MOFA <- create_mofa(GC_MOFA)
DataOptions <- get_default_data_options(GC_MOFA)
ModelOptions <- get_default_model_options(GC_MOFA)
ModelOptions$num_factors <- 5
TrainOptions <- get_default_training_options(GC_MOFA)
#TrainOptions$drop_factor_threshold <- 2019
#TrainOptions$DropFactorThreshold <- 0.01
#TrainOptions$verbose <- TRUE

GC_MOFA <- prepare_mofa(GC_MOFA, model_options = ModelOptions)

#Rogue_MOFA@TrainOptions$tolerance <-  0.1

#Rogue_MOFA@TrainOptions$DropFactorThreshold <- 0

#Rogue_MOFA@ModelOptions$numFactors <-  5

#use_python(/usr/bin/python3.7, required=TRUE)
MOFA_results <- run_mofa(GC_MOFA)
save(MOFA_results, file="~/analysis_tmp/Peters_workshop/MOFA_object.RData")

plot_variance_explained(MOFA_results)
library(scales)
library(plyr)
cols <- as.character(mapvalues(sobj@meta.data$Celltype, sort(unique(sobj@meta.data$Celltype)), hue_pal()(8)))

pairs(MOFA_results@expectations$Z$group1, pch=20, cex=0.8, col=cols, main="MOFA factors, data integration")

#plot(MOFA_results@expectations$Z$group1[,c(2,4)], pch=20, cex=0.8, col=cols, xlab="Factor 2", ylab="Factor 4")
legend("topright", unique(pd$coarseCellTypes), fill=hue_pal()(10), ncol=2)

#Genes most correlated with ATAC
plot_weights(object = MOFA_results, view = 'ATAC', factors = 1, nfeatures = 3, abs = FALSE)

#Peaks most correlated with EXPR
plot_weights(object = MOFA_results, view = 'EXPR', factors = 4, nfeatures = 3, abs = FALSE)



############################
#MOFA deconvolute which peaks are most predictive of gene expression?

#Factor 1
expr_weights <- get_weights(MOFA_results)$EXPR
atac_weights <- get_weights(MOFA_results)$ATAC

atac_ranges <- GRanges(sub("-", ":", rownames(atac_weights)))
values(atac_ranges) <- atac_weights

#Import gene coordinates from Ensembl
options(timeout=10000)
gtf <- as.data.frame(rtracklayer::import('https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz'))
gtf <- gtf[gtf$type=="gene",]
gtf$seqnames <- paste0("chr", gtf$seqnames)
gtfranges <- makeGRangesFromDataFrame(gtf)
values(gtfranges) <- gtf[,10:14]
names(gtfranges) <- gtf$gene_name
gtfranges <- gtfranges[names(gtfranges) %in% rownames(expr_weights)]


#Filter out nonoverlapping chroms
commonchroms <- intersect(seqnames(atac_ranges), seqnames(gtfranges))
gtfranges <- gtfranges[seqnames(gtfranges) %in% commonchroms]
atac_ranges <- atac_ranges[seqnames(atac_ranges) %in% commonchroms]



#Nearest gene to peak
atac_ranges$distToExprGene <- distanceToNearest(atac_ranges, gtfranges)
atac_ranges$nearestGene <- gtfranges$gene_name[subjectHits(atac_ranges$distToExprGene)]
atac_ranges$distToExprGene <- atac_ranges$distToExprGene@elementMetadata$distance

boxplot(atac_ranges$Factor1 ~ atac_ranges$distToExprGene < 10000, xlab="Is this peak within 10KB of an expressed gene?", ylab="Coefficient of Factor 1")

ol <- findOverlaps(atac_ranges + 10000, gtfranges)
res <- as.data.frame(ol)
res$atacpeak <- gsub(":", "-", paste(atac_ranges[queryHits(ol)]))
res$gene <- names(gtfranges[subjectHits(ol)])
res$Factor1dist <- sqrt((atac_weights[res$atacpeak, "Factor1"])^2 + (expr_weights[res$gene, "Factor1"])^2)
res$Factor2dist <- sqrt((atac_weights[res$atacpeak, "Factor2"])^2 + (expr_weights[res$gene, "Factor2"])^2)
res$Factor3dist <- sqrt((atac_weights[res$atacpeak, "Factor3"])^2 + (expr_weights[res$gene, "Factor3"])^2)
res$Factor4dist <- sqrt((atac_weights[res$atacpeak, "Factor4"])^2 + (expr_weights[res$gene, "Factor4"])^2)
res$Factor5dist <- sqrt((atac_weights[res$atacpeak, "Factor5"])^2 + (expr_weights[res$gene, "Factor5"])^2)

res <- res[order(res$Factor1dist, decreasing = T),]

resforplot <- data.frame(ATAC.coef.1=atac_weights[res$atacpeak, "Factor1"],
                   ATAC.coef.2=atac_weights[res$atacpeak, "Factor2"],
                   ATAC.coef.3=atac_weights[res$atacpeak, "Factor3"],
                   ATAC.coef.4=atac_weights[res$atacpeak, "Factor4"],
                   ATAC.coef.5=atac_weights[res$atacpeak, "Factor5"],
                   EXPR.coef.1=expr_weights[res$gene, "Factor1"],
                   EXPR.coef.2=expr_weights[res$gene, "Factor2"],
                   EXPR.coef.3=expr_weights[res$gene, "Factor3"],
                   EXPR.coef.4=expr_weights[res$gene, "Factor1"],
                   EXPR.coef.5=expr_weights[res$gene, "Factor5"],
                   atacpeak=res$atacpeak,
                   gene=res$gene,
                   F1dist=res$Factor1dist,
                   F2dist=res$Factor2dist,
                   F3dist=res$Factor3dist,
                   F4dist=res$Factor4dist,
                   F5dist=res$Factor5dist
                   )

resforplot$label <- ""
resforplot$label[1:10] <- resforplot$gene[1:10] 

library(ggrepel)
ggplot(resforplot) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_vline(xintercept=0, linetype="dashed", color = "black") + 
  geom_point(aes(x = ATAC.coef.1, y = EXPR.coef.1), colour="grey70") +
  geom_text_repel(aes(x = ATAC.coef.1, y = EXPR.coef.1, label = label),force = 2) +
  xlab("ATAC coefficients, Factor 1") + 
  ylab("RNA coefficients, Factor 1") +
  scale_color_manual(values=c('grey70')) +
  theme_linedraw() +
  annotate("label", x = -0.4, y = 0.5, label = "r=0.36") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1.25))) 

#Factor 3
resforplot <- resforplot[order(resforplot$F3dist, decreasing = T),]
resforplot$label <- ""
resforplot$label[1:10] <- resforplot$gene[1:10] 
library(ggrepel)
ggplot(resforplot) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_vline(xintercept=0, linetype="dashed", color = "black") + 
  geom_point(aes(x = ATAC.coef.3, y = EXPR.coef.3), colour="grey70") +
  geom_text_repel(aes(x = ATAC.coef.3, y = EXPR.coef.3, label = label),force = 2) +
  xlab("ATAC coefficients, Factor 3") + 
  ylab("RNA coefficients, Factor 3") +
  scale_color_manual(values=c('grey70')) +
  theme_linedraw() +
  annotate("label", x = -0.2, y = 0.5, label = "r=0.35") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1.25))) 

#Factor 4
resforplot <- resforplot[order(resforplot$F4dist, decreasing = T),]
resforplot$label <- ""
resforplot$label[1:10] <- resforplot$gene[1:10] 
library(ggrepel)
ggplot(resforplot) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  geom_vline(xintercept=0, linetype="dashed", color = "black") + 
  geom_point(aes(x = ATAC.coef.4, y = EXPR.coef.4), colour="grey70") +
  geom_text_repel(aes(x = ATAC.coef.4, y = EXPR.coef.4, label = label),force = 2) +
  xlab("ATAC coefficients, Factor 4") + 
  ylab("RNA coefficients, Factor 4") +
  scale_color_manual(values=c('grey70')) +
  theme_linedraw() +
  annotate("label", x = -0.2, y = 0.5, label = "r=0.35") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = rel(1.25))) 

#BCL11B
bcl11b <- res[res$gene=="BCL11B",]
bcl11batac <- logCPM_ATAC[bcl11b$atacpeak,]
bcl11batac <- rbind(bcl11batac, logCPM_RNA["BCL11B",])
rownames(bcl11batac)[6] <- "BCL11B_expr"
bcl11batac  <- as.matrix(bcl11batac)

typecol <- cols
names(typecol) <- colnames(bcl11batac)
anno.colors <- list(Celltype=c(`B cells`="#F8766D", `CD4+ T cells`="#CD9600",
                               `CD8+ T cells`="#7CAE00", `Dendritic cells`= "#00BE67",
                               Monocytes="#00BFC4", `NK cells`= "#00A9FF", Progenitors="#C77CFF", `T cells`="#FF61CC"))
anno <- data.frame(Celltype=sobj$Celltype)
rownames(anno) <- colnames(bcl11batac)
bcl11batac <- bcl11batac[c(6, 1:5),]
library(pheatmap)
pheatmap(bcl11batac, fontsize_col = 5, clustering_method = "ward.D", 
         annotation=anno, annotation_colors=anno.colors, show_colnames = F, cluster_rows = F)

#cd74
cd74 <- res[res$gene=="CD74",]
cd74atac <- logCPM_ATAC[cd74$atacpeak,]
cd74atac <- rbind(cd74atac, logCPM_RNA["CD74",])
rownames(cd74atac)[3] <- "cd74_expr"
cd74atac  <- as.matrix(cd74atac)

typecol <- cols
names(typecol) <- colnames(cd74atac)
anno.colors <- list(Celltype=c(`B cells`="#F8766D", `CD4+ T cells`="#CD9600",
                               `CD8+ T cells`="#7CAE00", `Dendritic cells`= "#00BE67",
                               Monocytes="#00BFC4", `NK cells`= "#00A9FF", Progenitors="#C77CFF", `T cells`="#FF61CC"))
anno <- data.frame(Celltype=sobj$Celltype)
rownames(anno) <- colnames(cd74atac)
cd74atac <- cd74atac[c(3, 1:2),]
library(pheatmap)
pheatmap(cd74atac, fontsize_col = 5, clustering_method = "ward.D", 
         annotation=anno, annotation_colors=anno.colors, show_colnames = F, cluster_rows = F)

