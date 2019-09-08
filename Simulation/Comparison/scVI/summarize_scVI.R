# Apply scVI to the simulated study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
library(ggplot2)

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

set.seed(12345)
# Working directory
setwd("/your/working/directory/BUSseq_implementation-1.0/Simulation/Comparison/scVI/")

###################
# Load simulation Data #
###################
# load the dimentsion information
dim <- unlist(read.table("../../RawCountData/dim_simulation_v4.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Loading the true cell type indicators
metadata <- read.table("../../RawCountData/metadata_simulation_v4.txt")

##############################
# load the inference by scVI #
##############################
w_scVI <- unlist(read.table("scVI_simulation_v4_clusters.txt"))
scVI_corrected <- read.table("scVI_simulation_v4_latent.txt")

#######
# ARI #
#######
ARI_scVI <- adjustedRandIndex(metadata$celltype, w_scVI)

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
scVI_dist <- dist(scVI_corrected) 
all.dists.scVI <- as.matrix(scVI_dist)
tsne_scVI_dist <- Rtsne(all.dists.scVI, is_distance=TRUE, perplexity = 30)

scVI_by_celltype<- "Image/tsne_simulation_scVI_by_celltype.jpeg"
plot_by_celltype(scVI_by_celltype, tsne_scVI_dist$Y)

scVI_by_batch <- "Image/tsne_simulation_scVI_by_batch.jpeg"
plot_by_batch(scVI_by_batch, tsne_scVI_dist$Y)

##################
# Draw PCA plots #
##################
pca.scVI <- prcomp(scVI_corrected, rank=2)
scVI_PCA<- "Image/pca_simulation_scVI.jpeg"
plot_by_celltype(scVI_PCA, pca.scVI$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
# save.image("scVI_workspace.RData")
save(ARI_scVI,tsne_scVI_dist,file = "scVI_results.RData")