# Apply scVI to the hematopoietic study.
rm(list=ls())
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot
set.seed(123)

# Working directory
setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/MouseHematopoietic/Comparison/scVI/")

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
w_scVI <- unlist(read.table("scVI_hemat_v1_clusters_0716.txt"))
scVI_corrected <- read.table("scVI_hemat_v1_latent_0716.txt")

#######
# ARI #
#######
ARI_scVI <- adjustedRandIndex(metadata$CellType, w_scVI)

#######################################################################
# Silhouette coefficient by estimated cell type labels of each method #
#######################################################################
# Silhouette coefficients of the latent variables
scVI_dist <- dist(scVI_corrected) 
sil_scVI <- silhouette(as.integer(w_scVI), dist = scVI_dist)
sil_scVI_true <- silhouette(as.integer(metadata$CellType), dist = scVI_dist)

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
all.dists.scVI <- as.matrix(scVI_dist)
tsne_scVI_dist <- Rtsne(all.dists.scVI, is_distance=TRUE, perplexity = 30)

scVI_by_celltype<- "Image/tsne_hemat_scVI_by_celltype.jpeg"
plot_by_celltype(scVI_by_celltype, tsne_scVI_dist$Y)

scVI_by_batch <- "Image/tsne_hemat_scVI_by_batch.jpeg"
plot_by_batch(scVI_by_batch, tsne_scVI_dist$Y)

##################
# Draw PCA plots #
##################
pca.scVI <- prcomp(scVI_corrected, rank=2)
scVI_PCA<- "Image/pca_hemat_scVI.jpeg"
plot_by_celltype(scVI_PCA, pca.scVI$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
save.image("scVI_workspace.RData")