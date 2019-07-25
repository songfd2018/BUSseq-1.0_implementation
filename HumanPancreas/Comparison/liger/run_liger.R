# library(devtools)
# install_github('MacoskoLab/liger')
rm(list=ls())
library(liger)
library(mclust) # For ARI
library(cluster) # For Silhouette
library(Rtsne) # For t-SNE plot

# coloring
library(scales)
require(WGCNA)
library(RColorBrewer)

set.seed(12345)
# Working directory
# setwd("/scratch/data01/BUSseq_cpp/BUSseq_implementation_v1/HumanPancreas/Comparison/liger/")
setwd("/your/working/directory/BUSseq_implementation-1.0/HumanPancreas/Comparison/liger/")

###################
# Load pancreas Data #
###################
# Loading pancreas count data
countdata <- read.table("../../RawCountData/count_data_pancreas_v1.txt") 

# Load dimension
dim <- unlist(read.table("../../RawCountData/dim_pancreas_v1.txt"))
N <- dim[1]
G <- dim[2]
B <- dim[3]
nb <- dim[3 + 1:B]

# Load metadata
metadata <- read.table("../../RawCountData/metadata_pancreas_v1.txt")

# Load gene_list
gene_list <- unlist(read.table("../../RawCountData/gene_list_pancreas_v1.txt",stringsAsFactors = FALSE))

rownames(countdata) <- gene_list
colnames(countdata) <- rownames(metadata)
pancreas.list <- list(GSE81076 = countdata[,1:nb[1]], GSE85241 = countdata[,1:nb[2]+nb[1]], GSE86473 = countdata[,1:nb[3]+sum(nb[1:2])], EMTAB5061 = countdata[,1:nb[4]+sum(nb[1:3])])

################################################################################################################
# Conduct analysis referring to https://github.com/MacoskoLab/liger/blob/master/vignettes/walkthrough_pbmc.pdf #
################################################################################################################
# Create liger object
liger_pancreas <- createLiger(pancreas.list, remove.missing = FALSE)

# Normalization
liger_pancreas <- normalize(liger_pancreas)
liger_pancreas@var.genes <- gene_list
liger_pancreas <- scaleNotCenter(liger_pancreas)

# running suggestK on multiple cores can greatly decrease the runtime
# looking for the section of the plot where this metric stops increasing as sharply
k.suggest <- suggestK(liger_pancreas, num.cores = 5, gen.new = T, return.data = T, plot.log2 = F)


# Take the lowest objective of three factorizations with different initializations
# Multiple restarts are recommended for initial analyses since iNMF is non-deterministic
liger_pancreas <- optimizeALS(liger_pancreas, k=30, thresh = 5e-5, nrep = 3)

liger_pancreas <- runTSNE(liger_pancreas, use.raw = T)
p1 <- plotByDatasetAndCluster(liger_pancreas, return.plots = T)
# Plot by dataset
print(p1[[1]])


##########################
# Currently working here #
##########################

# To better integrate the datasets, we perform a quantile alignment step. This process first identifies similarly
# loading cells across datasets by building a similarity graph based on shared factor neighborhoods.
liger_pancreas <- quantileAlignSNF(liger_pancreas, resolution = 0.4, small.clust.thresh = 20)

# liger_pancreas <- runTSNE(liger_pancreas)
# p_a <- plotByDatasetAndCluster(liger_pancreas, return.plots = T)
# # Modify plot output slightly
# p_a[[1]] <- p_a[[1]] + theme_classic() + theme(legend.position = c(0.85, 0.15)) +
#   guides(col=guide_legend(title = '', override.aes = list(size = 4)))
# print(p_a[[1]])
# clusters_prior <- readRDS('pbmc_alignment/tenx_seqwell_clusters.RDS')
# tenx_c <- droplevels(clusters_prior[rownames(a.pbmc@H[[1]])])
# seqwell_c <- droplevels(clusters_prior[rownames(a.pbmc@H[[2]])])
# metadata <- read.table(file = "../../RawCountData/pancreas_metadata.txt")
# GSE72857_celltype <- droplevels(metadata$CellType[1:nb[1]])
# GSE81682_celltype <- droplevels(metadata$CellType[1:nb[2] + nb[1]])
# # set_node_order = list(c(2, 6, 3, 4, 8, 1, 5, 7), c(1, 6, 2, 3, 4, 5, 7, 8, 9),
# #                       c(5, 2, 3, 6, 1, 4))
# # makeRiverplot(liger_pancreas, GSE72857_celltype, GSE81682_celltype, min.frac = 0.05,
# # node.order = set_node_order,
# # river.usr = c(0, 1, -0.6, 1.6)
# 
# print(p_a[[2]])
# 
# # This function uses the V loadings to identify dataset-specific genes and a
# # combination of the W and V loadings to identify shared markers
# # identify some shared and dataset-specific markers for each factor and plot them to help in cluster
# # annotation.
# markers <- getFactorMarkers(liger_pancreas, dataset1='GSE72857', dataset2='GSE81682',
#                             num.genes = 10)
# marker_genes <- markers$shared[order(markers$shared$p_value),]
w_liger <- liger_pancreas@alignment.clusters

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
loading_liger <- liger_pancreas@H.norm
liger_dist <- dist(loading_liger)
all.dists.liger <- as.matrix(liger_dist)
tsne_liger_dist <- Rtsne(all.dists.liger, is_distance=TRUE, perplexity = 30)

liger_by_celltype<- "Image/tsne_pancreas_liger_by_celltype.jpeg"
plot_by_celltype(liger_by_celltype, tsne_liger_dist$Y)

liger_by_batch <- "Image/tsne_pancreas_liger_by_batch.jpeg"
plot_by_batch(liger_by_batch, tsne_liger_dist$Y)

##################
# Draw PCA plots #
##################
pca.liger <- prcomp(loading_liger, rank=2)
liger_PCA<- "Image/pca_pancreas_liger.jpeg"
plot_by_celltype(liger_PCA, pca.liger$x, xlab = "PC 1", ylab = "PC 2")

# Store the workspace
#save.image("liger_workspace.RData")
save(ARI_liger,tsne_liger_dist,file = "liger_results.RData")