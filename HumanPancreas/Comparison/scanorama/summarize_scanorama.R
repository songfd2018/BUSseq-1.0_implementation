#Apply scanorama to the pancreatic study.
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
# setwd("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/HumanPancreas/Comparison/scanorama/")
setwd("/your/working/directory/BUSseq_implementation-1.0/HumanPancreas/Comparison/scanorama/")

###################
# Load pancreas Data #
###################
# load the dimentsion information
dim <- unlist(read.table("../../RawCountData/dim_pancreas_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Loading the true cell type indicators
metadata <- read.table("../../RawCountData/metadata_pancreas_v1.txt")

##############################
# load the inference by scVI #
##############################
scanorama_embedding <- read.table("scanorama_pancreas_v1_integrated.txt")
scanorama_corrected <- read.table("scanorama_pancreas_v1_corrected.txt")

clu_scanorama <- pam(scanorama_corrected, 7)
w_scanorama <- clu_scanorama$clustering

#######
# ARI #
#######
ARI_scanorama <- adjustedRandIndex(metadata$CellType, w_scanorama)

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
scanorama_dist <- dist(scanorama_embedding) 
all.dists.scanorama <- as.matrix(scanorama_dist)
tsne_scanorama_dist <- Rtsne(all.dists.scanorama, is_distance=TRUE, perplexity = 30)

scanorama_by_celltype<- "Image/tsne_pancreas_scanorama_by_celltype.jpeg"
plot_by_celltype(scanorama_by_celltype, tsne_scanorama_dist$Y)

scanorama_by_batch <- "Image/tsne_pancreas_scanorama_by_batch.jpeg"
plot_by_batch(scanorama_by_batch, tsne_scanorama_dist$Y)

##################
# Draw PCA plots #
##################
pca.scanorama <- prcomp(scanorama_embedding, rank=2)
scanorama_PCA<- "Image/pca_pancreas_scanorama.jpeg"
plot_by_celltype(scanorama_PCA, pca.scanorama$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
# save.image("scanorama_workspace.RData")
save(ARI_scanorama,tsne_scanorama_dist,file = "scanorama_results.RData")