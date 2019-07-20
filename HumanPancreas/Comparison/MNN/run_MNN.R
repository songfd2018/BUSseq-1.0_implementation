#Apply MNN to the pancreas study.
rm(list=ls())
library(scran)
library(cluster)

# coloring
require(WGCNA)
library(RColorBrewer)


set.seed(12345)

########################
# Load Pancreas Data #
########################
# Working directory
setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/HumanPancreas/Comparison/MNN")

# Loading the file name list of all pancreas count data
load("../../RawCountData/pancreas_countdata.RData")

B <- length(PancreasCounts)
nb <- rep(NA,B)
for(b in 1:B){
  nb[b] <- ncol(PancreasCounts[[b]])
}

# Load metadata
metadata <- read.table("../../RawCountData/pancreas_metadata.txt")

##################################
# Apply MNN to the Pancreas Data #
##################################
data_MNN <- PancreasCounts
# Referring to https://github.com/MarioniLab/MNN2017/blob/master/Pancreas
# Referring to http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html
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
  ave[[b]] <- max(1e-6,ave[[b]])
  if(b>1){
    data_MNN_normalized[[b]] <- data_MNN_normalized[[b]] * median(ave[[1]]/ave[[b]])
  }
}

# MNN batch correction
Xmnn <- mnnCorrect(data_MNN_normalized[[1]],
                   data_MNN_normalized[[2]],
                   data_MNN_normalized[[3]],
                   data_MNN_normalized[[4]],
                   svd.dim=0,
                   cos.norm.in=TRUE, cos.norm.out=TRUE,
                   var.adj=TRUE, 
                   k=20, sigma=0.1)

# combine corrected matrices together
corrected.df <- do.call(cbind.data.frame, Xmnn$corrected)
corrected.mat <- as.matrix(t(corrected.df))

##############
# Clustering #
##############

omat <- do.call(cbind, data_MNN_normalized)
sce <- SingleCellExperiment(list(logcounts=log1p(omat)))
reducedDim(sce, "Corrected") <- corrected.mat
sce$Batch <- rep(paste0("Batch",1:B),nb)

start_time <- Sys.time()
snn.gr <- buildSNNGraph(sce, use.dimred="Corrected")
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
MNN_dist <- dist(corrected.mat)
sil_MNN <- silhouette(as.integer(w_MNN), dist = MNN_dist)
sil_MNN_true <- silhouette(as.integer(metadata$CellType), dist = MNN_dist)

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
all.dists.MNN <- as.matrix(MNN_dist)
tsne_MNN_dist <- Rtsne(all.dists.MNN, is_distance=TRUE, perplexity = 30)

MNN_by_celltype<- "Image/tsne_pancreas_MNN_by_celltype.jpeg"
plot_by_celltype(MNN_by_celltype, tsne_MNN_dist$Y)

MNN_by_batch <- "Image/tsne_pancreas_MNN_by_batch.jpeg"
plot_by_batch(MNN_by_batch, tsne_MNN_dist$Y)

##########################################
# Draw the PCA plot on common cell types #
##########################################
pca.MNN <- prcomp(corrected.mat, rank=2)
MNN_PCA<- "Image/pca_pancreas_MNN.jpeg"
plot_by_celltype(MNN_PCA, pca.MNN$x, xlab = "PC 1", ylab = "PC 2")

# legend
pdf(file="Image/legend.pdf", width=10, height=8)
plot(c(-5,5),c(-4,4),type="n", bty="n", axes=FALSE, xlab="", ylab="")
legend(x=-5, y=4, legend=c("GSE81076","GSE85241","GSE86473","E-MTAB-5061"), pch=20, cex=2.5, col=alpha(colors4,0.6), title = "Batch", bty="n")
legend(x=-0, y=4, legend=names(table(metadata$CellType)), pch=20, cex=2.5, col=alpha(c("turquoise","blue","brown","orange1","green","deeppink","black"),0.6), title = "Cell Type", bty="n")
dev.off()

# Store the workspace
save.image("MNN_workspace.RData")

