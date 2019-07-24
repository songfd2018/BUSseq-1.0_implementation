# library(devtools)
# install_github('MacoskoLab/liger')
rm(list=ls())
library(liger)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot

# Working directory
#setwd("F:/scRNA/code/0601Cpp_BUSseq/BUSseq_implementation_v1/MouseHematopoietic/Comparison/liger/")
setwd("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/MouseHematopoietic/Comparison/liger/")

###################
# Load hemat Data #
###################
# Loading hematopoietic count data
countdata <- read.table("../../RawCountData/count_data_hemat_v1.txt") 

dim <- unlist(read.table("../../RawCountData/dim_hemat_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]
# Load metadata
metadata <- read.table("../../RawCountData/metadata_hemat_v1.txt")

# Load gene_list
gene_list <- unlist(read.table("../../RawCountData/gene_list_hemat_v1.txt",stringsAsFactors = FALSE))

rownames(countdata) <- gene_list
colnames(countdata) <- paste0("Batch",rep(1:2,nb),"_Cell",c(1:nb[1],1:nb[2]))
hemat.list <- list(GSE72857 = countdata[,1:nb[1]], GSE81682 = countdata[,1:nb[2]+nb[1]])


################################################################################################################
# Conduct analysis referring to https://github.com/MacoskoLab/liger/blob/master/vignettes/walkthrough_pbmc.pdf #
################################################################################################################
# Create liger object
liger_hemat <- createLiger(hemat.list, remove.missing = FALSE)

# Normalization
liger_hemat <- normalize(liger_hemat)
liger_hemat@var.genes <- gene_list
liger_hemat <- scaleNotCenter(liger_hemat)

# running suggestK on multiple cores can greatly decrease the runtime
# looking for the section of the plot where this metric stops increasing as sharply
# k.suggest <- suggestK(liger_hemat, num.cores = 5, gen.new = T, return.data = T, plot.log2 = F)

# Take the lowest objective of three factorizations with different initializations
# Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
liger_hemat <- optimizeALS(liger_hemat, k=22, thresh = 5e-5, nrep = 3)

liger_hemat <- runTSNE(liger_hemat, use.raw = T)
p1 <- plotByDatasetAndCluster(liger_hemat, return.plots = T)
# Plot by dataset
print(p1[[1]])

# To better integrate the datasets, we perform a quantile alignment step. This process first identifies similarly
# loading cells across datasets by building a similarity graph based on shared factor neighborhoods.
liger_hemat <- quantileAlignSNF(liger_hemat, resolution = 0.4, small.clust.thresh = 20)

# liger_hemat <- runTSNE(liger_hemat)
# p_a <- plotByDatasetAndCluster(liger_hemat, return.plots = T)
# # Modify plot output slightly
# p_a[[1]] <- p_a[[1]] + theme_classic() + theme(legend.position = c(0.85, 0.15)) +
#   guides(col=guide_legend(title = '', override.aes = list(size = 4)))
# print(p_a[[1]])

# clusters_prior <- readRDS('pbmc_alignment/tenx_seqwell_clusters.RDS')
# tenx_c <- droplevels(clusters_prior[rownames(a.pbmc@H[[1]])])
# seqwell_c <- droplevels(clusters_prior[rownames(a.pbmc@H[[2]])])


metadata <- read.table(file = "../../RawCountData/hemat_metadata.txt")
GSE72857_celltype <- droplevels(metadata$CellType[1:nb[1]])
GSE81682_celltype <- droplevels(metadata$CellType[1:nb[2] + nb[1]])
# set_node_order = list(c(2, 6, 3, 4, 8, 1, 5, 7), c(1, 6, 2, 3, 4, 5, 7, 8, 9),
#                       c(5, 2, 3, 6, 1, 4))
# makeRiverplot(liger_hemat, GSE72857_celltype, GSE81682_celltype, min.frac = 0.05,
              # node.order = set_node_order,
              # river.usr = c(0, 1, -0.6, 1.6)

print(p_a[[2]])

# This function uses the V loadings to identify dataset-specific genes and a
# combination of the W and V loadings to identify shared markers
# identify some shared and dataset-specific markers for each factor and plot them to help in cluster
# annotation.
markers <- getFactorMarkers(liger_hemat, dataset1='GSE72857', dataset2='GSE81682',
                            num.genes = 10)
marker_genes <- markers$shared[order(markers$shared$p_value),]
w_liger <- liger_hemat@alignment.clusters

#######
# ARI #
#######
ARI_liger <- adjustedRandIndex(metadata$CellType,w_liger)

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
loading_liger <- liger_hemat@H.norm
liger_dist <- dist(loading_liger)
all.dists.liger <- as.matrix(liger_dist)
tsne_liger_dist <- Rtsne(all.dists.liger, is_distance=TRUE, perplexity = 30)

liger_by_celltype<- "Image/tsne_hemat_liger_by_celltype.jpeg"
plot_by_celltype(liger_by_celltype, tsne_liger_dist$Y)

liger_by_batch <- "Image/tsne_hemat_liger_by_batch.jpeg"
plot_by_batch(liger_by_batch, tsne_liger_dist$Y)

##################
# Draw PCA plots #
##################
pca.liger <- prcomp(loading_liger, rank=2)
liger_PCA<- "Image/pca_hemat_liger.jpeg"
plot_by_celltype(liger_PCA, pca.liger$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
save.image("liger_workspace.RData")