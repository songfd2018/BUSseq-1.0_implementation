# Apply ZINB-WaVE to the pancreas study.
# Referring to https://github.com/drisso/zinb_analysis/blob/master/real_data/allen_covariates_1000.Rmd
rm(list=ls())
library(zinbwave)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
set.seed(12345)# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

########################
# Load Simulation Data #
########################
# Working directory
setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/HumanPancreas/Comparison/ZINBWaVE")

# Loading pancreas count data
load("../../RawCountData/pancreas_countdata.RData")

B <- length(PancreasCounts)
nb <- rep(NA, B)
for(b in 1:B){
  nb[b] <- ncol(PancreasCounts[[b]])
}

metadata <- read.table("../../RawCountData/pancreas_metadata.txt")


########################################
# Apply ZINB-WaVE to the Pancreas Data #
########################################

data_ZINBW <- NULL
for(b in 1:B){
  data_ZINBW <- cbind(data_ZINBW, PancreasCounts[[b]])
}

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

#######################################################################
# Silhouette coefficient by estimated cell type labels of each method #
#######################################################################
# Silhouette coefficients of the 10 components inferred by ZINBWaVE
ZINBWaVE_dist <- dist(zinb_batch@W)
sil_ZINBWaVE <- silhouette(as.integer(w_ZINBWaVE), dist = ZINBWaVE_dist)
sil_ZINBWaVE_true <-silhouette(as.integer(metadata$CellType), dist = ZINBWaVE_dist)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

allcolors<-labels2colors(metadata$CellType)#[allsamples])

allcolors[allcolors=="red"]<-"deeppink"
allcolors[allcolors=="yellow"]<-"orange1"#"darkgoldenrod1"

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=20, 
       col=alpha(allcolors[subset],0.6),
       xlab=xlab, ylab=ylab, main=main) 
  dev.off()
}

colors4<-brewer.pal(4,"Set3")
batch.cols<-colors4[rep(1:B,nb)]
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(6,6,4,2),cex.axis=2,cex.main=3,cex.lab=2.5)
  plot(Y[,1], Y[,2], cex=2,
       pch=20, 
       col=alpha(batch.cols[subset],0.6),
       xlab=xlab, ylab=ylab, main=main)#,  xlab="tSNE 1",ylab="tSNE 2")
  dev.off()
}

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
all.dists.ZINBW <- as.matrix(ZINBWaVE_dist)
tsne_ZINBW_dist <- Rtsne(all.dists.ZINBW, is_distance=TRUE, perplexity = 30)

ZINBW_by_celltype<- "Image/tsne_pancreas_ZINBWaVE_by_celltype.jpeg"
plot_by_celltype(ZINBW_by_celltype, tsne_ZINBW_dist$Y)

ZINBW_by_batch <- "Image/tsne_pancreas_ZINBWaVE_by_batch.jpeg"
plot_by_batch(ZINBW_by_batch, tsne_ZINBW_dist$Y)

##################
# Draw PCA plots #
##################
pca.ZINBW <- prcomp(zinb_batch@W, rank=2)
ZINBW_PCA<- "./Image/pca_pancreas_ZINBWaVE.jpeg"
plot_by_celltype(ZINBW_PCA, pca.ZINBW$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
save.image("ZINBWaVE_workspace.RData")