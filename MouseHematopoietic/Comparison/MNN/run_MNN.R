#Apply MNN to the hematopoietic study.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("scran")

rm(list=ls())
library(scran)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot

set.seed(123)

###################
# Load hemat Data #
###################
# Working directory
#setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/MouseHematopoietic/Comparison/MNN")
setwd("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/Comparison/MNN")

# Loading hematopoietic count data
load("../../RawCountData/hemat_countdata.RData")
HematCounts <- list(GSE72857 = dataA2,
                    GSE81682 = dataF2)
B <- 2
nb <- c(ncol(dataA2), ncol(dataF2))

# Load metadata
metadata <- read.table("../../RawCountData/hemat_metadata.txt")

#######################################
# Apply MNN to the Hematopoietic Data #
#######################################
data_MNN <- HematCounts
# Referring to https://github.com/MarioniLab/MNN2017/blob/master/Haematopoiesis/plotCorrections.R
# Normalization
data_MNN_normalized<-list()
for(b in 1:B){
  high.abF <- scater::calcAverage(data_MNN[[b]]) > 1
  clustF <- quickCluster(data_MNN[[b]], min.size=10 , method="igraph", subset.row=high.abF)
  sizeF <- computeSumFactors(data_MNN[[b]], sizes=seq(11, 81, 5), cluster=clustF, subset.row=high.abF)
  data_MNN_normalized[[b]] <- t(t(data_MNN[[b]])/sizeF)
}

# Rescaling the first dataset to match the coverage of the second.
ave<-list()
for(b in 1:B){
  ave[[b]] <- rowMeans(data_MNN_normalized[[b]])
  if(b>1){
    data_MNN_normalized[[b]] <- data_MNN_normalized[[b]] * median(ave[[1]]/ave[[b]])
  }
}

# Performing log-transformation and save results to file.
log_data_MNN<-list()
for(b in 1:B){
  log_data_MNN[[b]] <- log1p(data_MNN_normalized[[b]])
}

# Performing the correction with MNN (turned down the sigma to improve mixing).

mnn.out<-mnnCorrect(log_data_MNN[[1]], log_data_MNN[[2]],k=20, sigma=0.1,cos.norm.in=TRUE, cos.norm.out=TRUE, var.adj=TRUE,compute.angle=TRUE)

X.mnn <- cbind(mnn.out$corrected[[1]], mnn.out$corrected[[2]])
t.mnn <- t(X.mnn)

##############
# Clustering #
##############
# Referring to http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html
omat <- do.call(cbind, log_data_MNN)
sce <- SingleCellExperiment(list(logcounts=omat))
reducedDim(sce, "MNN") <- t.mnn
sce$Batch <- rep(paste0("Batch",1:2),nb)

start_time <- Sys.time()
snn.gr <- buildSNNGraph(sce, use.dimred="MNN")
clusters <- igraph::cluster_walktrap(snn.gr)
end_time <- Sys.time()
time_consumption <- end_time - start_time

table(clusters$membership, sce$Batch)

# The estimated cell type indicators by MNN
w_MNN <- factor(clusters$membership)

#######
# ARI #
#######
ARI_MNN <- adjustedRandIndex(metadata$CellType,w_MNN)

#######################################################################
# Silhouette coefficient by estimated cell type labels of each method #
#######################################################################
MNN_dist <- dist(t.mnn)
sil_MNN <- silhouette(as.integer(w_MNN), dist = MNN_dist)
sil_MNN_true <- silhouette(as.integer(metadata$CellType), dist = MNN_dist)

#################################
# scatter plots (t-SNE and PCA) #
#################################
if(!dir.exists("Image")){
  dir.create("Image")
}
# coloring for cell types
color.legendF <- c(MEP="orange", GMP="chartreuse4", CMP="magenta", 
                   LTHSC="dodgerblue", MPP="blue", LMPP="light blue", other="grey")
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
all.dists.mnn <- as.matrix(MNN_dist)
tsne_MNN_dist <- Rtsne(all.dists.mnn, is_distance=TRUE, perplexity = 30)

MNN_by_celltype<- "Image/tsne_hemat_MNN_by_celltype.jpeg"
plot_by_celltype(MNN_by_celltype, tsne_MNN_dist$Y)

MNN_by_batch <- "Image/tsne_hemat_MNN_by_batch.jpeg"
plot_by_batch(MNN_by_batch, tsne_MNN_dist$Y)

##################
# Draw PCA plots #
##################
pca.mnn <- prcomp(t.mnn, rank=2)
MNN_PCA<- "Image/pca_hemat_MNN.jpeg"
plot_by_celltype(MNN_PCA, pca.mnn$x, xlab = "PC 1", ylab = "PC 2")


# Store the workspace
save.image("MNN_workspace.RData")
