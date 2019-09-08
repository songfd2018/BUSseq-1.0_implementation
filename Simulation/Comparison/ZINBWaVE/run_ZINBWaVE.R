# Apply ZINB-WaVE to the simulation study.
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

# Working directory
setwd("/your/working/directory/BUSseq_implementation-1.0/Simulation/Comparison/ZINBWaVE/")

########################
# Load Simulated Data #
########################
# Loading the file name list of all simulation count data
countdata <- read.table("../../RawCountData/count_data_simulation_v4.txt") 

dim <- unlist(read.table("../../RawCountData/dim_simulation_v4.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_simulation_v4.txt")
gene_list <- paste0("gene_",1:G)

colnames(countdata) <- rownames(metadata)
rownames(countdata) <- gene_list


##########################################
# Apply ZINB-WaVE to the Simulatied Data #
##########################################
data_ZINBW <- as.matrix(countdata)

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
ARI_ZINBWaVE <- adjustedRandIndex(metadata$celltype,w_ZINBWaVE)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}

#####set cell type colorings
# batch color bar
color_by_batch<-c("#EB4334","#FBBD06","#35AA53","#4586F3")

# cell type color bar
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE")


# plot colored by cell type
allcolors<-color_by_celltype[metadata$celltype]#[allsamples])

plot_by_celltype<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  # The image with legend
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(10,10,2,2))
  plot(Y[,1], Y[,2],t="n",xaxt="n",yaxt="n",xlab="", ylab="")
  axis(1,line = 2.5, cex.axis=6,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  mtext(xlab, side=1, line=7.5, cex=6)
  mtext(ylab, side=2, line=5, cex=6)
  points(Y[,1], Y[,2], cex=3,
         pch=20, 
         col=alpha(allcolors[subset],0.6)) 
  dev.off()
}

# plot colored by batch
batch.cols<-color_by_batch[rep(1:B,nb)]
plot_by_batch<-function(pic_name,Y,subset=NULL,...,xlab = "tSNE 1", ylab = "tSNE 2",main=""){
  if (is.null(subset)) {
    subset <- seq_len(nrow(Y))
  }
  jpeg(pic_name,width = 1440, height = 1080, quality = 100)
  par(mfrow=c(1,1),mar=c(10,10,2,2))
  plot(Y[,1], Y[,2],t="n",xaxt="n",yaxt="n",xlab="", ylab="")
  axis(1,line = 2.5, cex.axis=6,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  mtext(xlab, side=1, line=7.5, cex=6)
  mtext(ylab, side=2, line=5, cex=6)
  points(Y[,1], Y[,2], cex=3,
         pch=20, 
         col=alpha(batch.cols[subset],0.6)) 
  dev.off()
}

#################################################
# Draw t-SNE plots by batch and true cell types #
#################################################
set.seed(123)
ZINBWaVE_dist <- dist(zinb_batch@W)
all.dists.ZINBW <- as.matrix(ZINBWaVE_dist)
tsne_ZINBW_dist <- Rtsne(all.dists.ZINBW, is_distance=TRUE, perplexity = 30)

ZINBW_by_celltype<- "Image/tsne_simulation_ZINBWaVE_by_celltype.jpeg"
plot_by_celltype(ZINBW_by_celltype, tsne_ZINBW_dist$Y)

ZINBW_by_batch <- "Image/tsne_simulation_ZINBWaVE_by_batch.jpeg"
plot_by_batch(ZINBW_by_batch, tsne_ZINBW_dist$Y)

##################
# Draw PCA plots #
##################
pca.ZINBW <- prcomp(zinb_batch@W, rank=2)
ZINBW_PCA<- "./Image/pca_simulation_ZINBWaVE.jpeg"
plot_by_celltype(ZINBW_PCA, pca.ZINBW$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
save.image("ZINBWaVE_workspace.RData")
save(ARI_ZINBWaVE,tsne_ZINBW_dist,file = "ZINBWaVE_results.RData")