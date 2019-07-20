#Apply scanorama to the hematopoietic study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot

# Working directory
setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/MouseHematopoietic/Comparison/scanorama/")

###################
# Load hemat Data #
###################
# load the dimentsion information
dim <- unlist(read.table("../../RawCountData/dim_hemat_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[1:B + 3]

# Loading the true cell type indicators
metadata <- read.table("../../RawCountData/hemat_metadata.txt")

##############################
# load the inference by scVI #
##############################
scanorama_embedding <- read.table("scanorama_hemat_v1_integrated.txt")
scanorama_corrected <- read.table("scanorama_hemat_v1_corrected.txt")

clu_scanorama <- pam(scanorama_corrected, 7)
w_scanorama <- clu_scanorama$clustering

#######
# ARI #
#######
ARI_scanorama <- adjustedRandIndex(metadata$CellType, w_scanorama)

#######################################################################
# Silhouette coefficient by estimated cell type labels of each method #
#######################################################################
# Silhouette coefficients of the low dimensional embedding space
scanorama_dist <- dist(scanorama_embedding) 
sil_scanorama <- silhouette(as.integer(w_scanorama), dist = scanorama_dist)
sil_scanorama_true <- silhouette(as.integer(metadata$CellType), dist = scanorama_dist)

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
all.dists.scanorama <- as.matrix(scanorama_dist)
tsne_scanorama_dist <- Rtsne(all.dists.scanorama, is_distance=TRUE, perplexity = 30)

scanorama_by_celltype<- "./Image/tsne_hemat_scanorama_by_celltype.jpeg"
plot_by_celltype(scanorama_by_celltype, tsne_scanorama_dist$Y)

scanorama_by_batch <- "./Image/tsne_hemat_scanorama_by_batch.jpeg"
plot_by_batch(scanorama_by_batch, tsne_scanorama_dist$Y)

##################
# Draw PCA plots #
##################
pca.scanorama <- prcomp(scanorama_embedding, rank=2)
scanorama_PCA<- "./Image/pca_hemat_scanorama.jpeg"
plot_by_celltype(scanorama_PCA, pca.scanorama$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
save.image("scanorama_workspace.RData")