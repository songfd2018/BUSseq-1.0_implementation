# Apply ZINB-WaVE to the hematopoietic study.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("zinbwave")

# Referring to https://github.com/drisso/zinb_analysis/blob/master/real_data/allen_covariates_1000.Rmd
rm(list=ls())
library(zinbwave)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

set.seed(123)

###########################
# Load Hematopoietic Data #
###########################
# Working directory
setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/MouseHematopoietic/Comparison/ZINBWaVE")
#setwd("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/Comparison/ZINBWaVE")

# Loading hematopoietic count data
data_ZINBW <- as.matrix(read.table("../../RawCountData/count_data_hemat_v1.txt"))

# Load dimension
dim <- read.table("../../RawCountData/dim_hemat_v1.txt")
dim <- unlist(dim)
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3+1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_hemat_v1.txt")
rownames(metadata) <- paste0("Batch_",rep(c(1,2),nb),"_Cell_",c(1:nb[1],1:nb[2]))
gene_list <- unlist(read.table("../../RawCountData/gene_list_hemat_v1.txt",stringsAsFactors = FALSE))

colnames(data_ZINBW) <- paste0("Batch_",rep(c(1,2),nb),"_Cell_",c(1:nb[1],1:nb[2]))
rownames(data_ZINBW) <- gene_list

#############################################
# Apply ZINB-WaVE to the Hematopoietic Data #
#############################################

# Factorizing the batch indicators
batch_ind <- factor(rep(1:B,nb))

# Performing ZINB-WaVE
zinb_batch <- zinbFit(data_ZINBW, K = 10, X=model.matrix(~batch_ind), epsilon=1e3)

# Clustering
# Referring to https://bioconductor.org/packages/release/bioc/vignettes/zinbwave/inst/doc/intro.html
# library("clusterExperiment")
# 
# ZINBW_clust<-clusterSingle(t(zinb_batch@W),sequential= T, subsample =F,mainClusterArgs=list(clusterFunction="kmeans"),seqArgs=list(k0=5,beta=0.95))
# w_ZINBW<-ZINBW_clust@clusterMatrix

data_ZINBW <- SummarizedExperiment(assays = data_ZINBW)
assayNames(data_ZINBW)[1] <- "counts"
merged_zinb <- zinbwave(data_ZINBW, fitted_model = zinb_batch, K = 10, epsilon=1000)

library(Seurat)

seu <- as.Seurat(x = merged_zinb, counts = "counts", data = "counts")

seu
seu <- FindNeighbors(seu, reduction = "zinbwave",
                     dims = 1:10 #this should match K
)
seu <- FindClusters(object = seu)

w_ZINBWaVE <- seu$seurat_clusters

#######
# ARI #
#######
ARI_ZINBWaVE <- adjustedRandIndex(metadata$CellType,w_ZINBWaVE)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}
# coloring for cell types
color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta", 
                   HSPC="cyan", LTHSC="dodgerblue", MPP="blue", LMPP="light blue", other="grey")
colmatF <- col2rgb(color.legendF) 
colmatA <- colmatF + 100 # A lighter shade.
colmatA[colmatA > 255] <- 255
color.legendA <- setNames(rgb(colmatA[1,], colmatA[2,], colmatA[3,], maxColorValue=255), names(color.legendF))
allcolors <- c(color.legendA[as.character(metadata$CellType[1:nb[1]])],color.legendF[as.character(metadata$CellType[1:nb[2]+nb[1]])])
first.batch <- rep(c(TRUE, FALSE), nb)

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 1, 21)[subset], 
       col=ifelse(first.batch, allcolors, "black")[subset],
       bg=allcolors[subset], xlab=xlab, ylab=ylab, main=main) 
  dev.off()
}

batchcolor <- c("lightcoral","lavender")
batch_ind <- rep(1:2,nb)
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=ifelse(first.batch, 1, 21)[subset], 
       col=ifelse(first.batch, batchcolor[batch_ind[subset]], "black"),
       bg=batchcolor[batch_ind[subset]], xlab=xlab, ylab=ylab, main=main)#,  xlab="tSNE 1",ylab="tSNE 2")
  dev.off()
}

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
ZINBWaVE_dist <- dist(zinb_batch@W)
all.dists.ZINBW <- as.matrix(ZINBWaVE_dist)
tsne_ZINBW_dist <- Rtsne(all.dists.ZINBW, is_distance=TRUE, perplexity = 30)

ZINBW_by_celltype<- "Image/tsne_hemat_ZINBWaVE_by_celltype.jpeg"
plot_by_celltype(ZINBW_by_celltype, tsne_ZINBW_dist$Y)

ZINBW_by_batch <- "Image/tsne_hemat_ZINBWaVE_by_batch.jpeg"
plot_by_batch(ZINBW_by_batch, tsne_ZINBW_dist$Y)

##################
# Draw PCA plots #
##################
pca.ZINBW <- prcomp(zinb_batch@W, rank=2)
ZINBW_PCA<- "./Image/pca_hemat_ZINBWaVE.jpeg"
plot_by_celltype(ZINBW_PCA, pca.ZINBW$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
save.image("ZINBWaVE_workspace.RData")